"""Rebuild daily SQLite examples from msPurityData's regression-v1 directory.

Run from the msPurity checkout with the regression directory as the argument.
Keep original identifiers, both runs, linked scans, and a scan from each file
even when a feature has no links. The matching tests also exercise ties across
the full set of averaged query spectra.
"""
from pathlib import Path
import shutil
import sqlite3
import sys
import tempfile


def compact(source, target):
    with tempfile.TemporaryDirectory(dir=target.parent) as tmp:
        output = Path(tmp) / target.name
        shutil.copyfile(source, output)
        db = sqlite3.connect(output)
        try:
            tables = {r[0] for r in db.execute(
                "SELECT name FROM sqlite_master WHERE type='table'")}
            groups = {12, 27, 53, 89, 162, 187, 410}
            if source.name == "createDatabase_example.sqlite":
                groups.update(int(r[0]) for r in db.execute(
                    "SELECT DISTINCT grpid FROM s_peak_meta WHERE grpid IS NOT NULL"))
            if "xcms_match" in tables:
                groups.update(r[0] for r in db.execute("SELECT grpid FROM xcms_match"))
            db.execute("CREATE TEMP TABLE keep_groups (id INTEGER PRIMARY KEY)")
            db.executemany("INSERT INTO keep_groups VALUES (?)", [(x,) for x in groups])
            db.executescript("""
                DELETE FROM c_peak_X_c_peak_group WHERE grpid NOT IN keep_groups;
                DELETE FROM c_peak_groups WHERE grpid NOT IN keep_groups;
                DELETE FROM c_peaks WHERE cid NOT IN (SELECT cid FROM c_peak_X_c_peak_group);
                DELETE FROM c_peak_X_s_peak_meta WHERE cid NOT IN (SELECT cid FROM c_peaks);
                CREATE TEMP TABLE keep_meta (id INTEGER PRIMARY KEY);
                INSERT INTO keep_meta SELECT DISTINCT pid FROM c_peak_X_s_peak_meta;
                INSERT OR IGNORE INTO keep_meta SELECT MIN(pid) FROM s_peak_meta
                    WHERE fileid IS NOT NULL AND seqNum IS NOT NULL GROUP BY fileid;
            """)
            if "av_peaks" in tables:
                db.execute("DELETE FROM av_peaks WHERE grpid NOT IN keep_groups")
            else:
                db.execute("INSERT OR IGNORE INTO keep_meta SELECT pid FROM s_peak_meta "
                           "WHERE CAST(grpid AS INTEGER) IN keep_groups")
            # Historical scan peaks can have an incorrect stored pid. Keep the
            # metadata for both identifiers, so conversion can repair ownership.
            db.execute("DELETE FROM s_peaks WHERE pid NOT IN keep_meta")
            if "av_peaks" not in tables:
                db.execute("DELETE FROM s_peaks WHERE CAST(grpid AS INTEGER) NOT IN keep_groups")
            db.execute("INSERT OR IGNORE INTO keep_meta SELECT m.pid FROM s_peak_meta m "
                       "JOIN s_peaks p ON m.fileid = p.fileid AND m.seqNum = p.scan")
            db.execute("DELETE FROM s_peak_meta WHERE pid NOT IN keep_meta")
            if "av_peaks" not in tables:
                assert not db.execute(
                    "SELECT pid FROM s_peak_meta WHERE CAST(grpid AS INTEGER) "
                    "NOT IN (SELECT grpid FROM c_peak_groups)").fetchone()
            db.commit()
            assert db.execute("PRAGMA integrity_check").fetchone() == ("ok",)
            assert not db.execute("PRAGMA foreign_key_check").fetchall()
            db.execute("VACUUM")
        finally:
            db.close()
        output.replace(target)


if __name__ == "__main__":
    root = Path(sys.argv[1])
    for name in ("db/createDatabase_example.sqlite", "db/createDatabase_example_OLD.sqlite",
                 "db/create_database_example.sqlite", "sm/spectralMatching_result.sqlite"):
        compact(root / name, Path("inst/extdata/tests") / name)
