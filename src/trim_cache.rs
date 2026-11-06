use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

use anyhow::{Context, Result};

use crate::barcode::{cdna_trim_bounds, PutativeRow};

#[derive(Clone, Copy, Debug)]
pub(crate) struct CachedTrimBounds {
    pub start_5p: usize,
    pub end_3p: usize,
}

impl CachedTrimBounds {
    pub fn from_putative(row: &PutativeRow, read_len: usize) -> Self {
        Self {
            start_5p: cdna_trim_bounds(row, read_len, false, true)
                .map(|(start, _)| start)
                .unwrap_or(read_len),
            end_3p: cdna_trim_bounds(row, read_len, true, false)
                .map(|(_, end)| end)
                .unwrap_or(0),
        }
    }

    pub fn bounds(self, read_len: usize, valid_3p: bool, valid_5p: bool) -> Option<(usize, usize)> {
        let start = if valid_5p { self.start_5p } else { 0 };
        let end = if valid_3p { self.end_3p } else { read_len };
        (start < end && end <= read_len).then_some((start, end))
    }

    pub fn read_from(reader: &mut (impl Read + ?Sized)) -> Result<Self> {
        let mut bytes = [0u8; 16];
        reader
            .read_exact(&mut bytes)
            .context("read cDNA trim cache")?;
        Ok(Self {
            start_5p: usize::try_from(u64::from_le_bytes(bytes[..8].try_into()?))?,
            end_3p: usize::try_from(u64::from_le_bytes(bytes[8..].try_into()?))?,
        })
    }
}

pub(crate) struct TrimCache {
    writer: BufWriter<File>,
}

impl TrimCache {
    pub fn new(directory: Option<&Path>) -> Result<Self> {
        let file = if let Some(directory) = directory {
            fs::create_dir_all(directory)?;
            tempfile::tempfile_in(directory)
        } else {
            tempfile::tempfile()
        }
        .context("create cDNA trim cache; set --ont-temp-dir or TMPDIR to writable scratch")?;
        Ok(Self {
            writer: BufWriter::with_capacity(1024 * 1024, file),
        })
    }

    pub fn write_entry(&mut self, entry: CachedTrimBounds) -> Result<()> {
        self.writer
            .write_all(&(entry.start_5p as u64).to_le_bytes())?;
        self.writer
            .write_all(&(entry.end_3p as u64).to_le_bytes())?;
        Ok(())
    }

    pub fn reader(&mut self) -> Result<BufReader<&mut File>> {
        self.writer.flush()?;
        let file = self.writer.get_mut();
        file.seek(SeekFrom::Start(0))?;
        Ok(BufReader::with_capacity(1024 * 1024, file))
    }
}
