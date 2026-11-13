#!/usr/bin/env python3
"""
Compare the old crossing file with the new one made from the corrected database.

Author: Saurabh Shukla
Email: sshukla2@alaska.edu

Crossings are matched by date, flight and station. For the matched ones it
reports how much the crossing time, the closest distance and the speed
changed, and it lists the flights whose old speed was 0 with their new speed.
Crossings that exist in only one of the two files are counted too (a
corrected spike can move a flight closer to, or further from, a station).

Writes crossing_changes.csv (one row per matched crossing) and
zero_speed_fixed.csv next to the new crossing file. Nothing is changed.
Run from the parkshwynodal folder:  python compare_crossings.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

# ===== CONFIG =====
OLD = Path("input/all_station_crossing_db_UTM_20km_corrected_accepted_spikes.txt")
NEW = Path("input/all_station_crossing_db_UTM_20km_corrected_db.txt")
START_DATE = "20190211"   # only compare these dates (use the same dates as the new file)
END_DATE = "20190326"
# =====

COLS = ["Date", "FlightNum", "X_m", "Y_m", "HorizontalDistance_m", "ClosestTime", "Altitude_m",
        "Speed_mps", "Heading", "Station", "Equipment", "SlantDistance_km"]
KEY = ["Date", "FlightNum", "Station"]


def read(path):
    c = pd.read_csv(path, names=COLS, header=None, dtype={"Date": str, "FlightNum": str, "Station": str})
    c = c[(c.Date >= START_DATE) & (c.Date <= END_DATE)]
    return c.drop_duplicates(KEY)


def main():
    old, new = read(OLD), read(NEW)
    m = old.merge(new, on=KEY, how="outer", suffixes=("_old", "_new"), indicator=True)
    both = m[m["_merge"] == "both"].copy()

    both["dt_s"] = both.ClosestTime_new - both.ClosestTime_old
    both["dd_m"] = both.HorizontalDistance_m_new - both.HorizontalDistance_m_old
    both["dv_mps"] = both.Speed_mps_new - both.Speed_mps_old
    keep = KEY + ["ClosestTime_old", "ClosestTime_new", "dt_s", "HorizontalDistance_m_old",
                  "HorizontalDistance_m_new", "dd_m", "Speed_mps_old", "Speed_mps_new", "dv_mps"]
    both[keep].to_csv(NEW.parent / "crossing_changes.csv", index=False)

    zero = both[both.Speed_mps_old.fillna(0) == 0]
    zero[KEY + ["Speed_mps_old", "Speed_mps_new"]].to_csv(NEW.parent / "zero_speed_fixed.csv", index=False)

    def pct(x, q):
        return np.nanpercentile(np.abs(x), q) if len(x) else np.nan

    print(f"Dates {START_DATE} to {END_DATE}")
    print(f"Old crossings: {len(old)}   new crossings: {len(new)}")
    print(f"Matched: {len(both)}   only in old: {(m['_merge'] == 'left_only').sum()}   "
          f"only in new: {(m['_merge'] == 'right_only').sum()}")
    print("Change for matched crossings (median / 95th percentile / max of the absolute change):")
    for name, col, unit in (("crossing time", "dt_s", "s"), ("closest distance", "dd_m", "m"),
                            ("speed", "dv_mps", "m/s")):
        x = both[col]
        print(f"  {name:17s} {pct(x, 50):9.2f} / {pct(x, 95):9.2f} / {np.nanmax(np.abs(x)) if len(x) else np.nan:9.2f} {unit}")
    print(f"Crossings that changed time by more than 5 s: {(both.dt_s.abs() > 5).sum()}")
    print(f"Crossings with old speed 0: {len(zero)} from {zero.FlightNum.nunique()} flights; "
          f"new speed: median {zero.Speed_mps_new.median():.1f} m/s, "
          f"still 0: {(zero.Speed_mps_new.fillna(0) == 0).sum()}")
    print(f"Wrote {NEW.parent / 'crossing_changes.csv'} and {NEW.parent / 'zero_speed_fixed.csv'}")


if __name__ == "__main__":
    main()
