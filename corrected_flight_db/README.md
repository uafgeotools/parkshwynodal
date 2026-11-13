# Station crossings from the corrected flight path database

Author: Saurabh Shukla (sshukla2@alaska.edu)

| Script | What it does |
|---|---|
| `make_crossings_corrected.py` | Makes the station crossing file (12 columns). Default: the original FR24 files and reported speed, as before. With `--use-corrected`: the corrected flight path database (positions corrected, smoothed speed). |
| `compare_crossings.py` | Compares an old and a new crossing file: changes in crossing time, closest distance and speed, and the crossings that had speed 0. |

Run both from the `parkshwynodal` folder (they read `input/` and `src/` there):

    python corrected_flight_db/make_crossings_corrected.py --use-corrected --radius-km 20 --start-date 20190211 --end-date 20190326
    python corrected_flight_db/compare_crossings.py

The corrected database itself is made with the scripts in `corrected_flight_db/` of `parkshwynodal_supp`.
