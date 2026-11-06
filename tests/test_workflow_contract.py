import copy
import json
import unittest
from pathlib import Path


FIXTURE_DIR = Path(__file__).parent / "fixtures" / "cli"
COMMANDS = ("single", "mixed", "glycine", "analyze")
OPTION_REQUIRED_KEYS = {
    "name",
    "aliases",
    "type",
    "required",
    "required_when",
    "default",
    "conflicts",
    "description",
}
OPTION_ALLOWED_KEYS = OPTION_REQUIRED_KEYS | {"choices"}
POSITIONAL_KEYS = {
    "name",
    "type",
    "cardinality",
    "required",
    "description",
}
CONDITION_KEYS = {"option", "operator", "value"}
ACCEPTED_INVOCATION_KEYS = {"argv", "reason"}
REJECTED_INVOCATION_KEYS = {"argv", "reason", "violation"}
VIOLATION_KEYS = {"code", "option"}
ALLOWED_VIOLATION_CODES = {
    "command_mismatch",
    "conflict",
    "invalid_type",
    "invalid_value",
    "missing_conditional",
    "missing_positional",
    "missing_required",
    "missing_value",
    "unexpected_positional",
    "unexpected_value",
    "unknown_option",
}
ALLOWED_TYPES = {"bool", "float", "int", "path", "string", "string_list"}
CARDINALITIES = {"one", "zero_or_one", "one_or_more", "zero_or_more"}
COMMON_TOP_LEVEL_KEYS = {
    "schema_version",
    "id",
    "command",
    "description",
    "sources",
    "positionals",
    "options",
    "conditional_requirements",
    "execution",
    "invocations",
}
EXPECTED_TOP_LEVEL_KEYS = {
    "single": COMMON_TOP_LEVEL_KEYS
    | {"resolution_rules", "interpreter", "compatibility_exceptions"},
    "mixed": COMMON_TOP_LEVEL_KEYS
    | {
        "resolution_rules",
        "interpreter",
        "compatibility_exceptions",
        "singlet_threshold_semantics",
    },
    "glycine": COMMON_TOP_LEVEL_KEYS,
    "analyze": COMMON_TOP_LEVEL_KEYS | {"command_aliases", "resolution_rules"},
}


class ContractError(AssertionError):
    pass


def reject_duplicate_keys(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ContractError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def loads_contract(text):
    try:
        return json.loads(text, object_pairs_hook=reject_duplicate_keys)
    except json.JSONDecodeError as error:
        raise ContractError(f"invalid JSON: {error.msg}") from error


def load_contract(name):
    with (FIXTURE_DIR / f"{name}.json").open(encoding="utf-8") as handle:
        return loads_contract(handle.read())


def option_map(contract):
    return {option["name"]: option for option in contract["options"]}


def require_exact_keys(value, expected, label):
    if not isinstance(value, dict):
        raise ContractError(f"{label} must be an object")
    actual = set(value)
    if actual != expected:
        missing = sorted(expected - actual)
        extra = sorted(actual - expected)
        raise ContractError(f"{label} keys differ: missing={missing}, extra={extra}")


def require_string(value, label, *, nonempty=True):
    if not isinstance(value, str) or (nonempty and not value):
        raise ContractError(f"{label} must be a non-empty string")


def require_string_list(value, label, *, nonempty=False, allow_empty_items=False):
    if not isinstance(value, list) or (nonempty and not value):
        raise ContractError(f"{label} must be a list")
    for index, item in enumerate(value):
        if not isinstance(item, str) or (not allow_empty_items and not item):
            qualifier = "string" if allow_empty_items else "non-empty string"
            raise ContractError(f"{label}[{index}] must be a {qualifier}")


def value_matches_type(value, value_type):
    if value is None:
        return True
    if value_type == "bool":
        return type(value) is bool
    if value_type == "int":
        return type(value) is int
    if value_type == "float":
        return type(value) in (int, float)
    if value_type in ("path", "string"):
        return isinstance(value, str)
    if value_type == "string_list":
        return isinstance(value, list) and all(isinstance(item, str) for item in value)
    return False


def validate_option(option):
    if not isinstance(option, dict):
        raise ContractError("option must be an object")
    actual_keys = set(option)
    if not OPTION_REQUIRED_KEYS.issubset(actual_keys) or not actual_keys.issubset(
        OPTION_ALLOWED_KEYS
    ):
        missing = sorted(OPTION_REQUIRED_KEYS - actual_keys)
        extra = sorted(actual_keys - OPTION_ALLOWED_KEYS)
        raise ContractError(f"option keys differ: missing={missing}, extra={extra}")
    require_string(option["name"], "option.name")
    if not option["name"].startswith("--"):
        raise ContractError("option.name must start with --")
    require_string_list(option["aliases"], "option.aliases")
    if any(len(alias) < 2 or not alias.startswith("-") for alias in option["aliases"]):
        raise ContractError(f"{option['name']} aliases must start with -")
    if len(set(option["aliases"])) != len(option["aliases"]):
        raise ContractError(f"{option['name']} aliases must be unique")
    if option["type"] not in ALLOWED_TYPES:
        raise ContractError(f"unknown option type: {option['type']!r}")
    if type(option["required"]) is not bool:
        raise ContractError("option.required must be bool")
    if not isinstance(option["required_when"], list):
        raise ContractError("option.required_when must be a list")
    for condition in option["required_when"]:
        require_exact_keys(condition, CONDITION_KEYS, "condition")
        require_string(condition["option"], "condition.option")
        if condition["operator"] not in ("equals", "not_equals"):
            raise ContractError("condition.operator is invalid")
    require_string_list(option["conflicts"], "option.conflicts")
    if len(set(option["conflicts"])) != len(option["conflicts"]):
        raise ContractError(f"{option['name']} conflicts must be unique")
    require_string(option["description"], "option.description")
    if not value_matches_type(option["default"], option["type"]):
        raise ContractError(
            f"{option['name']} default is incompatible with type {option['type']}"
        )
    if option["required"] and option["default"] is not None:
        raise ContractError(f"{option['name']} required option cannot have a default")
    if "choices" in option:
        if not isinstance(option["choices"], list) or not option["choices"]:
            raise ContractError(f"{option['name']} choices must be a non-empty list")
        if any(
            not value_matches_type(choice, option["type"])
            for choice in option["choices"]
        ):
            raise ContractError(f"{option['name']} choice has incompatible type")
        if any(
            choice == earlier
            for index, choice in enumerate(option["choices"])
            for earlier in option["choices"][:index]
        ):
            raise ContractError(f"{option['name']} choices must be unique")
        if option["default"] is not None and option["default"] not in option["choices"]:
            raise ContractError(f"{option['name']} default is not an allowed choice")


def validate_positional(positional):
    require_exact_keys(positional, POSITIONAL_KEYS, "positional")
    require_string(positional["name"], "positional.name")
    if positional["name"].startswith("-"):
        raise ContractError("positional.name cannot start with -")
    if positional["type"] not in ALLOWED_TYPES - {"bool"}:
        raise ContractError("positional.type is invalid")
    if positional["cardinality"] not in CARDINALITIES:
        raise ContractError("positional.cardinality is invalid")
    if type(positional["required"]) is not bool:
        raise ContractError("positional.required must be bool")
    cardinality_is_required = positional["cardinality"] in ("one", "one_or_more")
    if positional["required"] != cardinality_is_required:
        raise ContractError("positional.required disagrees with cardinality")
    require_string(positional["description"], "positional.description")


def validate_invocation_schema(invocation, outcome):
    expected = (
        ACCEPTED_INVOCATION_KEYS if outcome == "accepted" else REJECTED_INVOCATION_KEYS
    )
    require_exact_keys(invocation, expected, f"{outcome} invocation")
    require_string_list(invocation["argv"], "invocation.argv", nonempty=True)
    require_string(invocation["reason"], "invocation.reason")
    if outcome == "rejected":
        require_exact_keys(invocation["violation"], VIOLATION_KEYS, "violation")
        require_string(invocation["violation"]["code"], "violation.code")
        require_string(invocation["violation"]["option"], "violation.option")
        if invocation["violation"]["code"] not in ALLOWED_VIOLATION_CODES:
            raise ContractError("violation.code is invalid")


def validate_top_level(contract, fixture_name):
    require_exact_keys(contract, EXPECTED_TOP_LEVEL_KEYS[fixture_name], "contract")
    if type(contract["schema_version"]) is not int or contract["schema_version"] != 1:
        raise ContractError("schema_version must be integer 1")
    if contract["id"] != fixture_name:
        raise ContractError("contract id does not match fixture name")
    require_string_list(contract["command"], "command", nonempty=True)
    if any(token.startswith("-") for token in contract["command"]):
        raise ContractError("command tokens cannot be options")
    require_string(contract["description"], "description")
    require_string_list(contract["sources"], "sources", nonempty=True)
    if not isinstance(contract["positionals"], list):
        raise ContractError("positionals must be a list")
    if not isinstance(contract["options"], list) or not contract["options"]:
        raise ContractError("options must be a non-empty list")
    require_string_list(contract["conditional_requirements"], "conditional_requirements")
    require_exact_keys(contract["execution"], {"stages", "does_not_require"}, "execution")
    require_string_list(contract["execution"]["stages"], "execution.stages")
    require_string_list(
        contract["execution"]["does_not_require"], "execution.does_not_require"
    )
    require_exact_keys(contract["invocations"], {"accepted", "rejected"}, "invocations")
    for outcome in ("accepted", "rejected"):
        if not isinstance(contract["invocations"][outcome], list) or not contract[
            "invocations"
        ][outcome]:
            raise ContractError(f"invocations.{outcome} must be a non-empty list")

    if "command_aliases" in contract:
        if not isinstance(contract["command_aliases"], list):
            raise ContractError("command_aliases must be a list")
        for alias in contract["command_aliases"]:
            require_string_list(alias, "command_alias", nonempty=True)
        aliases = [tuple(alias) for alias in contract["command_aliases"]]
        if len(set(aliases)) != len(aliases) or tuple(contract["command"]) in aliases:
            raise ContractError("command_aliases must be unique and differ from command")
    if "resolution_rules" in contract:
        require_string_list(contract["resolution_rules"], "resolution_rules")
    if "interpreter" in contract:
        interpreter = contract["interpreter"]
        require_exact_keys(
            interpreter, {"precedence", "environment", "fallback", "scope"}, "interpreter"
        )
        require_string_list(interpreter["precedence"], "interpreter.precedence", nonempty=True)
        for key in ("environment", "fallback", "scope"):
            require_string(interpreter[key], f"interpreter.{key}")
    if "compatibility_exceptions" in contract:
        compatibility = contract["compatibility_exceptions"]
        require_exact_keys(compatibility, {"rejected_options"}, "compatibility_exceptions")
        rejected_options = compatibility["rejected_options"]
        if not isinstance(rejected_options, dict) or not rejected_options:
            raise ContractError("rejected_options must be a non-empty object")
        for name, reason in rejected_options.items():
            require_string(name, "rejected option name")
            require_string(reason, "rejected option reason")
    if "singlet_threshold_semantics" in contract:
        semantics = contract["singlet_threshold_semantics"]
        require_exact_keys(
            semantics,
            {
                "parser",
                "accepted_values",
                "rejected_values",
                "range_constraints",
                "propagates_to",
            },
            "singlet_threshold_semantics",
        )
        require_string(semantics["parser"], "singlet parser")
        for key in ("accepted_values", "rejected_values"):
            require_string_list(
                semantics[key],
                f"singlet_threshold_semantics.{key}",
                allow_empty_items=True,
            )
        for key in ("range_constraints", "propagates_to"):
            require_string_list(semantics[key], f"singlet_threshold_semantics.{key}")


def validate_cross_references(contract):
    options = option_map(contract)
    if len(options) != len(contract["options"]):
        raise ContractError("duplicate canonical option name")
    spellings = {}
    positional_names = [positional["name"] for positional in contract["positionals"]]
    if len(set(positional_names)) != len(positional_names):
        raise ContractError("duplicate positional name")
    for index, positional in enumerate(contract["positionals"][:-1]):
        if positional["cardinality"] in ("one_or_more", "zero_or_more"):
            raise ContractError(
                f"variadic positional {positional['name']} must be last"
            )
    for option in contract["options"]:
        for spelling in (option["name"], *option["aliases"]):
            if spelling in spellings:
                raise ContractError(f"option spelling collision: {spelling}")
            spellings[spelling] = option["name"]
    for option in contract["options"]:
        for conflict in option["conflicts"]:
            if conflict not in options:
                raise ContractError(f"unknown conflict target: {conflict}")
            if option["name"] not in options[conflict]["conflicts"]:
                raise ContractError(f"asymmetric conflict: {option['name']} / {conflict}")
        for condition in option["required_when"]:
            if condition["option"] not in spellings:
                raise ContractError(f"unknown condition option: {condition['option']}")
            condition_name = spellings[condition["option"]]
            condition_option = options[condition_name]
            if not value_matches_type(condition["value"], condition_option["type"]):
                raise ContractError(
                    f"condition value for {condition['option']} has incompatible type"
                )


def parse_cli_value(raw_value, option):
    value_type = option["type"]
    try:
        if value_type == "int":
            return int(raw_value)
        if value_type == "float":
            return float(raw_value)
        if value_type == "string_list":
            return raw_value.split(",")
        if value_type in ("path", "string") and raw_value:
            return raw_value
    except ValueError:
        return None
    return None


def invocation_violations(contract, invocation):
    argv = invocation["argv"]
    command = contract["command"]
    if argv[: len(command)] != command:
        return [{"code": "command_mismatch", "option": "<command>"}]

    options = option_map(contract)
    spellings = {
        spelling: option["name"]
        for option in contract["options"]
        for spelling in (option["name"], *option["aliases"])
    }
    present = set()
    values = {name: option["default"] for name, option in options.items()}
    positionals = []
    violations = []
    tokens = argv[len(command) :]
    index = 0
    while index < len(tokens):
        token = tokens[index]
        if token.startswith("-"):
            spelling, separator, attached = token.partition("=")
            canonical = spellings.get(spelling)
            if canonical is None:
                violations.append({"code": "unknown_option", "option": spelling})
                index += 1
                continue
            option = options[canonical]
            present.add(canonical)
            if option["type"] == "bool":
                if separator:
                    violations.append({"code": "unexpected_value", "option": canonical})
                values[canonical] = True
                index += 1
                continue
            if separator:
                raw_value = attached
            elif index + 1 < len(tokens):
                raw_value = tokens[index + 1]
                index += 1
            else:
                violations.append({"code": "missing_value", "option": canonical})
                index += 1
                continue
            if option["type"] == "string_list" and not separator:
                raw_values = [raw_value]
                while index + 1 < len(tokens) and not tokens[index + 1].startswith("-"):
                    raw_values.append(tokens[index + 1])
                    index += 1
                parsed = raw_values
            else:
                parsed = parse_cli_value(raw_value, option)
            if parsed is None:
                violations.append({"code": "invalid_type", "option": canonical})
            elif "choices" in option and parsed not in option["choices"]:
                violations.append({"code": "invalid_value", "option": canonical})
            else:
                values[canonical] = parsed
            index += 1
            continue
        positionals.append(token)
        index += 1

    for option in contract["options"]:
        name = option["name"]
        if option["required"] and name not in present:
            violations.append({"code": "missing_required", "option": name})
        for condition in option["required_when"]:
            condition_name = spellings[condition["option"]]
            condition_matches = values[condition_name] == condition["value"]
            if condition["operator"] == "not_equals":
                condition_matches = not condition_matches
            if condition_matches and values[name] is None:
                violations.append({"code": "missing_conditional", "option": name})
        for conflict in option["conflicts"]:
            if name in present and conflict in present and name < conflict:
                violations.append({"code": "conflict", "option": name})

    positional_index = 0
    for positional in contract["positionals"]:
        remaining = len(positionals) - positional_index
        cardinality = positional["cardinality"]
        minimum = 1 if cardinality in ("one", "one_or_more") else 0
        maximum = 1 if cardinality in ("one", "zero_or_one") else remaining
        take = min(remaining, maximum)
        if take < minimum:
            violations.append(
                {"code": "missing_positional", "option": positional["name"]}
            )
        for raw_value in positionals[positional_index : positional_index + take]:
            parsed = parse_cli_value(raw_value, {"type": positional["type"]})
            if parsed is None:
                violations.append(
                    {"code": "invalid_type", "option": positional["name"]}
                )
        positional_index += take
    if positional_index < len(positionals):
        violations.append({"code": "unexpected_positional", "option": "<positionals>"})
    return violations


def validate_contract(contract, fixture_name):
    validate_top_level(contract, fixture_name)
    for option in contract["options"]:
        validate_option(option)
    for positional in contract["positionals"]:
        validate_positional(positional)
    validate_cross_references(contract)
    for outcome in ("accepted", "rejected"):
        for invocation in contract["invocations"][outcome]:
            validate_invocation_schema(invocation, outcome)
            violations = invocation_violations(contract, invocation)
            if outcome == "accepted" and violations:
                raise ContractError(f"accepted invocation has {violations[0]['code']}")
            if outcome == "rejected" and invocation["violation"] not in violations:
                raise ContractError(
                    f"rejected invocation does not exhibit {invocation['violation']}"
                )


class WorkflowContractTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contracts = {name: load_contract(name) for name in COMMANDS}

    def test_fixture_schema_is_complete(self):
        for fixture_name, contract in self.contracts.items():
            with self.subTest(fixture=fixture_name):
                validate_top_level(contract, fixture_name)
            for index, option in enumerate(contract["options"]):
                with self.subTest(
                    fixture=fixture_name,
                    option=option.get("name", "<missing>"),
                    index=index,
                ):
                    validate_option(option)
            for index, positional in enumerate(contract["positionals"]):
                with self.subTest(
                    fixture=fixture_name,
                    positional=positional.get("name", "<missing>"),
                    index=index,
                ):
                    validate_positional(positional)
            with self.subTest(fixture=fixture_name, section="cross_references"):
                validate_cross_references(contract)
            for outcome in ("accepted", "rejected"):
                for index, invocation in enumerate(contract["invocations"][outcome]):
                    with self.subTest(
                        fixture=fixture_name, outcome=outcome, index=index
                    ):
                        validate_invocation_schema(invocation, outcome)

    def test_option_names_and_aliases_are_unique_per_command(self):
        for fixture_name, contract in self.contracts.items():
            seen = {}
            for option in contract["options"]:
                for spelling in (option["name"], *option["aliases"]):
                    with self.subTest(fixture=fixture_name, spelling=spelling):
                        self.assertNotIn(spelling, seen)
                    seen[spelling] = option["name"]

    def test_aliases_defaults_conditions_and_conflicts_reference_known_options(self):
        for fixture_name, contract in self.contracts.items():
            options = option_map(contract)
            all_spellings = {
                spelling
                for option in contract["options"]
                for spelling in (option["name"], *option["aliases"])
            }
            for option in contract["options"]:
                with self.subTest(fixture=fixture_name, option=option["name"]):
                    for alias in option["aliases"]:
                        self.assertTrue(alias.startswith("-"))
                    for conflict in option["conflicts"]:
                        self.assertIn(conflict, options)
                        self.assertIn(option["name"], options[conflict]["conflicts"])
                    for condition in option["required_when"]:
                        self.assertIn(condition["option"], all_spellings)
                        self.assertIn(condition["operator"], {"equals", "not_equals"})
                        self.assertIn("value", condition)
                    if option["required"]:
                        self.assertIsNone(option["default"])

    def test_legacy_aliases_and_representative_defaults_are_frozen(self):
        workflow_aliases = {
            "--tso-seq": ["--tso_seq"],
            "--rtp-seq": ["--rtp_seq"],
            "--barcode-list-10bp": ["--barcode_list_10bp"],
            "--ref-dir": ["--ref_dir"],
            "--gene-fasta": ["--genome-fa"],
            "--out-dir": ["--outdir"],
            "--sample-id": ["--sample"],
            "--threads": ["--thread"],
            "--barcode-extract-mode": ["--barcode_extract_mode"],
            "--glycine-err": ["--err"],
            "--glycine-shift": ["--shift"],
            "--min-len": ["--min_len"],
            "--umi-len": ["--umi_len"],
            "--help": ["-h"],
            "--version": ["-V"],
        }
        for fixture_name in ("single", "mixed"):
            options = option_map(self.contracts[fixture_name])
            aliases = {
                name: option["aliases"]
                for name, option in options.items()
                if option["aliases"]
            }
            with self.subTest(fixture=fixture_name):
                self.assertEqual(workflow_aliases, aliases)
                self.assertEqual(32, options["--threads"]["default"])
                self.assertEqual(16, options["--cluster-threads"]["default"])
                self.assertEqual(10, options["--glycine-jobs"]["default"])
                self.assertEqual(64, options["--glycine-threads"]["default"])
                self.assertEqual("string_list", options["--fastq"]["type"])
                self.assertTrue(options["--light-output"]["default"])
                self.assertFalse(options["--full-output"]["default"])
                self.assertIsNone(options["--python"]["default"])

        glycine = option_map(self.contracts["glycine"])
        self.assertEqual("0.25,0.25", glycine["--err"]["default"])
        self.assertEqual("100,100", glycine["--shift"]["default"])
        self.assertEqual(100, glycine["--min_len"]["default"])
        self.assertEqual(10, glycine["--jobs"]["default"])
        self.assertEqual(64, glycine["--total-threads"]["default"])
        self.assertEqual("string_list", glycine["--fastq"]["type"])

        analyze = option_map(self.contracts["analyze"])
        self.assertEqual("Flora", analyze["--out_dir"]["default"])
        self.assertEqual(32, analyze["--threads"]["default"])
        self.assertEqual("fixed_seq", analyze["--barcode_extract_mode"]["default"])
        self.assertEqual(20, analyze["--min_reads_per_cell"]["default"])
        self.assertTrue(analyze["--absorb_unassigned_paired"]["default"])

    def test_legacy_toggle_pairs_are_ordered_overrides_not_conflicts(self):
        expected_rules = [
            "full-pipeline and upstream-only are accepted together; the last occurrence wins",
            "full-output and light-output are accepted together; the last occurrence wins",
        ]
        for fixture_name in ("single", "mixed"):
            contract = self.contracts[fixture_name]
            with self.subTest(fixture=fixture_name):
                self.assertEqual(expected_rules, contract["resolution_rules"][:2])
                conflicts = {
                    option["name"]: option["conflicts"]
                    for option in contract["options"]
                    if option["conflicts"]
                }
                self.assertEqual(
                    {"--fastq": ["--fastq-dir"], "--fastq-dir": ["--fastq"]},
                    conflicts,
                )

    def test_each_contract_has_representative_acceptance_and_rejection(self):
        for fixture_name, contract in self.contracts.items():
            for outcome in ("accepted", "rejected"):
                invocations = contract["invocations"][outcome]
                with self.subTest(fixture=fixture_name, outcome=outcome):
                    self.assertGreaterEqual(len(invocations), 2)
                for index, invocation in enumerate(invocations):
                    with self.subTest(
                        fixture=fixture_name, outcome=outcome, index=index
                    ):
                        violations = invocation_violations(contract, invocation)
                        if outcome == "accepted":
                            self.assertEqual([], violations)
                        else:
                            self.assertIn(invocation["violation"], violations)

    def test_full_workflows_freeze_python_precedence(self):
        for fixture_name in ("single", "mixed"):
            contract = self.contracts[fixture_name]
            options = option_map(contract)
            with self.subTest(fixture=fixture_name):
                self.assertIn("--python", options)
                self.assertEqual("path", options["--python"]["type"])
                self.assertEqual(
                    ["--python", "FLORA_PYTHON", "python3"],
                    contract["interpreter"]["precedence"],
                )
                self.assertEqual("FLORA_PYTHON", contract["interpreter"]["environment"])
                self.assertEqual("python3", contract["interpreter"]["fallback"])
                self.assertEqual(
                    "all_python_stages",
                    contract["interpreter"]["scope"],
                )

    def test_glycine_bin_dir_is_an_explicit_full_workflow_rejection(self):
        for fixture_name in ("single", "mixed"):
            contract = self.contracts[fixture_name]
            rejected = contract["compatibility_exceptions"]["rejected_options"]
            invocations = contract["invocations"]["rejected"]
            with self.subTest(fixture=fixture_name):
                self.assertEqual(
                    "must_error_never_ignore", rejected["--glycine-bin-dir"]
                )
                self.assertNotIn("--glycine-bin-dir", option_map(contract))
                self.assertTrue(
                    any("--glycine-bin-dir" in item["argv"] for item in invocations)
                )

    def test_mixed_singlet_threshold_preserves_python_float_semantics(self):
        contract = self.contracts["mixed"]
        option = option_map(contract)["--singlet-threshold"]
        semantics = contract["singlet_threshold_semantics"]
        self.assertEqual("float", option["type"])
        self.assertEqual(0.9, option["default"])
        self.assertEqual("python_float", semantics["parser"])
        self.assertEqual(
            ["barnyard_qc", "parameters_tsv"], semantics["propagates_to"]
        )
        self.assertEqual([], semantics["range_constraints"])

        for value in semantics["accepted_values"]:
            with self.subTest(value=value):
                float(value)
        self.assertEqual(-0.1, float("-0.1"))
        self.assertNotEqual(float("nan"), float("nan"))
        self.assertEqual(float("inf"), float("1e309"))
        self.assertEqual(float("-inf"), float("-1e309"))
        for value in semantics["rejected_values"]:
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    float(value)

        single = self.contracts["single"]
        self.assertNotIn("--singlet-threshold", option_map(single))

    def test_command_boundaries_are_isolated(self):
        glycine = self.contracts["glycine"]["execution"]
        self.assertEqual(["glycine"], glycine["stages"])
        self.assertEqual(
            ["reference", "alignment", "python"], glycine["does_not_require"]
        )

        analyze = self.contracts["analyze"]["execution"]
        self.assertEqual(
            ["barcode_assignment", "umi_assignment", "cell_assignment"],
            analyze["stages"],
        )
        self.assertEqual(
            ["glycine", "reference", "alignment", "python"],
            analyze["does_not_require"],
        )

        self.assertNotIn("--python", option_map(self.contracts["glycine"]))
        self.assertNotIn("--python", option_map(self.contracts["analyze"]))

    def test_unknown_option_mutation_is_rejected(self):
        contract = copy.deepcopy(self.contracts["single"])
        contract["invocations"]["accepted"][0]["argv"].append("--unknown")
        with self.assertRaisesRegex(ContractError, "unknown_option"):
            validate_contract(contract, "single")

    def test_bool_as_int_default_mutation_is_rejected(self):
        contract = copy.deepcopy(self.contracts["single"])
        option_map(contract)["--light-output"]["default"] = 1
        with self.assertRaisesRegex(ContractError, "default"):
            validate_contract(contract, "single")

    def test_changed_default_that_invalidates_an_example_is_rejected(self):
        contract = copy.deepcopy(self.contracts["single"])
        option_map(contract)["--skip-glycine"]["default"] = True
        with self.assertRaisesRegex(ContractError, "missing_conditional"):
            validate_contract(contract, "single")

    def test_positional_type_mutation_is_rejected(self):
        contract = copy.deepcopy(self.contracts["analyze"])
        contract["positionals"][0]["type"] = "int"
        with self.assertRaisesRegex(ContractError, "invalid_type"):
            validate_contract(contract, "analyze")

    def test_missing_option_value_mutation_is_rejected(self):
        contract = copy.deepcopy(self.contracts["single"])
        contract["invocations"]["accepted"][1]["argv"].pop()
        with self.assertRaisesRegex(ContractError, "missing_value"):
            validate_contract(contract, "single")

    def test_unknown_top_level_key_mutation_is_rejected(self):
        contract = copy.deepcopy(self.contracts["single"])
        contract["unexpected"] = True
        with self.assertRaisesRegex(ContractError, "extra=.*unexpected"):
            validate_contract(contract, "single")

    def test_unknown_option_key_mutation_is_rejected(self):
        contract = copy.deepcopy(self.contracts["single"])
        contract["options"][0]["unexpected"] = True
        with self.assertRaisesRegex(ContractError, "extra=.*unexpected"):
            validate_contract(contract, "single")

    def test_duplicate_json_key_is_rejected(self):
        with self.assertRaisesRegex(ContractError, "duplicate JSON key: id"):
            loads_contract('{"id": "single", "id": "mixed"}')


class WorkflowTimingContractTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[1]
        cls.scripts = {
            name: (root / name).read_text(encoding="utf-8")
            for name in ("run_all.sh", "run_all_mixed_species.sh")
        }

    def test_internal_stages_report_implementation_and_elapsed_time(self):
        for name, script in self.scripts.items():
            with self.subTest(script=name):
                self.assertIn('stage_impl="rust"', script)
                self.assertIn('stage_impl="python"', script)
            self.assertIn(
                'Stage ${rust_name} (${stage_impl}) elapsed:', script
            )
            self.assertIn('return "${stage_status}"', script)

    def test_metrics_workbook_and_compact_cleanup_are_wired_into_both_workflows(self):
        project_root = Path(__file__).resolve().parents[1]
        for runner in ("run_all.sh", "run_all_mixed_species.sh"):
            script = (project_root / runner).read_text(encoding="utf-8")
            self.assertIn(
                'run_stage metrics_summary "${DOWNSTREAM_DIR}/metrics_summary.py"',
                script,
            )
            self.assertIn('--output-xlsx "metrics_summary.xlsx"', script)
            self.assertIn(
                '--barcode-to-cell "${UPSTREAM_DIR}/barcode_to_cell.csv"', script
            )
            self.assertIn("cleanup_large_intermediates", script)
            self.assertIn(
                'if [[ "${LIGHT_OUTPUT}" -eq 1 && "${SAVE_INTERMEDIATE}" -eq 0 ]]',
                script,
            )
            self.assertIn("${SAMPLE_ID}.filtered.tagged.sorted.bam", script)

            cleanup = script.split("cleanup_large_intermediates() {", 1)[1].split(
                "\n}", 1
            )[0]
            for retained in (
                "filtered.tagged.sorted.bam",
                "gene_expression.tsv",
                "isoform_expression.tsv",
                "rna_cluster.tsv",
                "barcode_to_cell.csv",
                "metrics_summary.xlsx",
                "single_cell_report.html",
                "barnyard_qc",
            ):
                self.assertNotIn(retained, cleanup)

    def test_light_output_streams_alignment_fastq_without_persisting_it(self):
        for name, script in self.scripts.items():
            with self.subTest(script=name):
                self.assertIn("ALIGNMENT_STREAMING=0", script)
                self.assertIn(
                    'if [[ "${LIGHT_OUTPUT}" -eq 1 && "${SAVE_INTERMEDIATE}" -eq 0 ]]',
                    script,
                )
                self.assertIn("ALIGNMENT_STREAMING=1", script)
                self.assertIn('${SAMPLE_ID}.cdna.fastq.pipe', script)
                self.assertIn('mkfifo "${ALIGNMENT_FASTQ}"', script)
                self.assertIn('exec 9<>"${ALIGNMENT_FASTQ}"', script)
                self.assertIn('wait "${ALIGNMENT_STREAM_PID}"', script)
                self.assertIn('rm -f "${ALIGNMENT_FASTQ}"', script)
                self.assertIn('${SAMPLE_ID}.cdna.fastq.gz', script)

    def test_correction_maps_are_opt_in_debug_outputs(self):
        project_root = Path(__file__).resolve().parents[1]
        pipeline = (project_root / "src" / "pipeline.rs").read_text(encoding="utf-8")
        main = (project_root / "src" / "main.rs").read_text(encoding="utf-8")
        self.assertIn(
            "config.save_intermediate || config.save_correction_maps", pipeline
        )
        self.assertIn("barcode_correction.write_maps: skipped", pipeline)
        self.assertIn("save_correction_maps: cli.save_merge_debug", main)

    def test_remove_final_bam_is_opt_in_and_runs_only_after_the_report(self):
        for name, script in self.scripts.items():
            with self.subTest(script=name):
                self.assertIn("REMOVE_FINAL_BAM=0", script)
                self.assertIn(
                    "--remove-final-bam) REMOVE_FINAL_BAM=1; shift ;;", script
                )
                self.assertIn("cleanup_all_bam_outputs()", script)
                all_bam_cleanup = script.split("cleanup_all_bam_outputs() {", 1)[
                    1
                ].split("\n}", 1)[0]
                for bam in (
                    "${ALIGN_DIR}/${SAMPLE_ID}.aligned.sorted.bam",
                    "${MATRIX_DIR}/${SAMPLE_ID}.filtered.cb_ur.gn.sorted.bam",
                    "${MATRIX_DIR}/${SAMPLE_ID}.filtered.tagged.sorted.bam",
                ):
                    self.assertIn(bam, all_bam_cleanup)
                self.assertLess(
                    script.index("run_stage build_report"),
                    script.index('if [[ "${REMOVE_FINAL_BAM}" -eq 1 ]]'),
                )

    def test_direct_expensive_stages_have_named_timing(self):
        expected = (
            "Stage glycine (rust) elapsed:",
            "Stage alignment_pipeline (external) elapsed:",
            "Stage aligned_bam_index (external) elapsed:",
            "Stage bamtobed (external) elapsed:",
            "Stage gene_bam_index (external) elapsed:",
            "Stage rna_cluster_analysis (python) elapsed:",
        )
        for name, script in self.scripts.items():
            for marker in expected:
                with self.subTest(script=name, marker=marker):
                    self.assertIn(marker, script)


class SaturationStageContractTest(unittest.TestCase):
    def test_saturation_uses_embedded_rust_compute_and_python_plot_only(self):
        root = Path(__file__).resolve().parents[1]
        runtime = (root / "src" / "workflow_runtime.rs").read_text(encoding="utf-8")
        self.assertIn('"saturation",', runtime)
        for runner_name in ["run_all.sh", "run_all_mixed_species.sh"]:
            runner = (root / runner_name).read_text(encoding="utf-8")
            self.assertIn(
                'run_stage saturation "${DOWNSTREAM_DIR}/Saturation.py"', runner
            )
            self.assertIn('--plot-existing-tsv "${SAMPLE_ID}.saturation.tsv"', runner)


if __name__ == "__main__":
    unittest.main()
