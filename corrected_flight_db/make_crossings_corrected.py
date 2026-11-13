#!/usr/bin/env python3
"""
Make the station crossing file, from the original FR24 files (default,
USE_CORRECTED = False) or from the corrected flight path database
(USE_CORRECTED = True, or --use-corrected).

Author: Saurabh Shukla
Email: sshukla2@alaska.edu

Based on the parkshwynodal crossing script; the functions it uses from
src/doppler_funcs.py are unchanged.

Without --use-corrected the script reads the raw FR24 files and the reported
speed, as the crossing script always did. With --use-corrected it reads the
corrected flight path database (correct_flight_paths.py) instead. Everything
else (closest point, crossing time, averaging, station time check, 12 output
columns) is the same in both modes.

What changes in the input with --use-corrected:
    positions   latitude_corrected, longitude_corrected (spikes, altitude
                errors and stray points in gaps are corrected)
    altitude    altitude_corrected_ft (still feet, avg_return converts to m)
    speed       speed_smoothed_knots (still knots, avg_return converts to m/s).
                This is the smoothed FR24 speed; where FR24 reported 0 while
                the aircraft was moving, it is the speed of the corrected path.
    rows        positions flagged endpoint_removed and rows without a position
                are left out; flights too short to check (fewer than three
                positions) are used as they are, same as before

Output (same 12 columns as before, no header):
    corrected mode: input/all_station_crossing_db_UTM_{radius}km_corrected_db.txt
    original mode:  input/all_station_crossing_db_UTM_{radius}km_original.txt
    and next to it *_path_speed.csv: crossings whose speed came from the
    corrected positions (FR24 reported 0, or the FR24 speed was replaced)

Run from the parkshwynodal folder (it reads input/ and src/ there), for example:
    python corrected_flight_db/make_crossings_corrected.py --use-corrected --radius-km 20 --start-date 20190211 --end-date 20190326
"""

import os
import sys
import argparse
import pyproj
import concurrent.futures
import pandas as pd
import numpy as np
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# the repository folder (the one with src/), also when this script sits in a subfolder
repo_root = Path(__file__).resolve().parent
while not (repo_root / "src").is_dir() and repo_root != repo_root.parent:
    repo_root = repo_root.parent

if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))

from src.doppler_funcs import (
    load_flight_file,
    time_check,
    find_closest_point,
    get_equip,
    closest_time_calc,
    avg_return,
    flight_lat_lon_to_utm,
)


NODE_FILE = Path("input/parkshwy_nodes.txt")
# USE_CORRECTED = False: the original FR24 files and the reported speed (the crossing file as before)
# USE_CORRECTED = True:  the corrected flight path database and the smoothed speed
# (the command line options --use-corrected / --use-original override this default)
USE_CORRECTED = False
RAW_DIR = Path("/scratch/irseppi/nodal_data/flightradar24")
CORRECTED_DB = Path.home() / "REPOSITORIES/specufex/ten_flight_velocity_study/corrected_flight_db_v1/flights"

# rows of the corrected database that should not be used for crossings
# (flights too short to check are kept as they are, like in the old crossing file;
# rows without any position are dropped below)
SKIP_STATUS = {"endpoint_removed"}

UTM_ZONE = "6"
utm_proj = pyproj.Proj(proj="utm", zone=UTM_ZONE, ellps="WGS84")


CROSSING_COLUMNS = [
    "Date",
    "FlightNum",
    "X_m",
    "Y_m",
    "HorizontalDistance_m",
    "ClosestTime",
    "Altitude_m",
    "Speed_mps",
    "Heading",
    "Station",
    "Equipment",
    "SlantDistance_km",
]


def parse_args():
    """
    Read the command line options.

    Returns:
        argparse.Namespace: Radius, workers, database folder, mode and dates.
    """
    parser = argparse.ArgumentParser(
        description="Generate all-aircraft node crossing file from the corrected flight path database."
    )

    parser.add_argument(
        "--radius-km",
        type=float,
        default=20.0,
        help="Maximum horizontal closest-approach distance in km. Example: 20 or 30.",
    )

    parser.add_argument(
        "--workers",
        type=int,
        default=os.cpu_count(),
        help="Number of parallel workers.",
    )

    parser.add_argument(
        "--corrected-db",
        type=Path,
        default=CORRECTED_DB,
        help="Folder with the corrected flight files (YYYYMMDD_flightid.csv).",
    )

    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--use-corrected",
        dest="use_corrected",
        action="store_true",
        default=USE_CORRECTED,
        help="Use the corrected flight path database and the smoothed speed.",
    )
    mode.add_argument(
        "--use-original",
        dest="use_corrected",
        action="store_false",
        help="Use the original FR24 files and the reported speed (as before).",
    )

    parser.add_argument(
        "--start-date",
        default="20190211",
        help="First date to include, YYYYMMDD.",
    )

    parser.add_argument(
        "--end-date",
        default="20190326",
        help="Last date to include, YYYYMMDD.",
    )

    return parser.parse_args()


def radius_label(radius_km):
    """
    Radius as text for file names (20.0 gives 20, 2.5 gives 2p5).

    Args:
        radius_km (float): Radius in km.

    Returns:
        str: The label.
    """
    if float(radius_km).is_integer():
        return str(int(radius_km))

    return str(radius_km).replace(".", "p")


def corrected_flight_list(corrected_db, start_date, end_date):
    """
    List the corrected flight files between two dates (both included).

    Args:
        corrected_db (pathlib.Path): Folder of the corrected database.
        start_date (str): First date, YYYYMMDD.
        end_date (str): Last date, YYYYMMDD.

    Returns:
        tuple: Paths of the flight files and their file names (two lists).
    """

    flight_files = []
    filenames = []

    for f in sorted(Path(corrected_db).glob("*.csv")):
        date = f.name[0:8]
        if date.isdigit() and start_date <= date <= end_date:
            flight_files.append(str(f))
            filenames.append(f.name)

    return flight_files, filenames


def original_flight_list(raw_dir, start_date, end_date):
    """
    List the original FR24 files between two dates (both included).

    Files filed in the folder of another day are skipped.

    Args:
        raw_dir (pathlib.Path): Folder of the original FR24 database.
        start_date (str): First date, YYYYMMDD.
        end_date (str): Last date, YYYYMMDD.

    Returns:
        tuple: Paths of the flight files and their file names (two lists).
    """

    flight_files = []
    filenames = []

    for folder in sorted(Path(raw_dir).glob("*_positions")):
        for f in sorted(folder.glob("*.csv")):
            date = f.name[0:8]
            if date.isdigit() and start_date <= date <= end_date and folder.name.startswith(date):
                flight_files.append(str(f))
                filenames.append(f.name)

    return flight_files, filenames


def load_original_flight_file(flight_file, filename):
    """
    Load one original FR24 file with the original loader.

    Args:
        flight_file (str): Path of the file.
        filename (str): File name, YYYYMMDD_flightid.csv.

    Returns:
        tuple: The values of load_flight_file and None for the speed source, or
        None when the flight has fewer than two positions.
    """

    loaded = load_flight_file(flight_file, filename)
    if len(loaded[3]) < 2:
        return None

    return (*loaded, None)


def load_corrected_flight_file(flight_file, filename):
    """
    Load one flight from the corrected database, in the units load_flight_file returns.

    Altitude stays in feet and speed in knots, because avg_return converts both.
    Flagged end points and rows without a position are dropped and the index is
    reset, because avg_return and closest_time_calc use positions index and
    index + 1.

    Args:
        flight_file (str): Path of the corrected file.
        filename (str): File name, YYYYMMDD_flightid.csv.

    Returns:
        tuple: UTM x and y (km), flight path, timestamps, altitude (ft), speed
        (knots), heading, flight number, date and speed source; or None when the
        flight has fewer than two positions.
    """

    flight_data = pd.read_csv(flight_file, sep=",", low_memory=False)

    flight_data = flight_data[~flight_data["position_status"].isin(SKIP_STATUS)]
    flight_data = flight_data.dropna(subset=["latitude_corrected", "longitude_corrected", "snapshot_id"])
    flight_data = flight_data.drop_duplicates("snapshot_id").reset_index(drop=True)

    # a path needs at least two positions; return None and the flight is skipped quietly
    if len(flight_data) < 2:
        return None

    flight_latitudes = flight_data["latitude_corrected"]
    flight_longitudes = flight_data["longitude_corrected"]
    timestamp = flight_data["snapshot_id"]
    alt = flight_data["altitude_corrected_ft"]
    speed = flight_data["speed_smoothed_knots"]
    head = flight_data["heading"]
    speed_source = flight_data["speed_source"]

    fname = filename
    flight_num = fname[9:18]
    date = fname[0:8]

    flight_utm_x_km, flight_utm_y_km, flight_path = flight_lat_lon_to_utm(flight_latitudes, flight_longitudes)

    return flight_utm_x_km, flight_utm_y_km, flight_path, timestamp, alt, speed, head, flight_num, date, speed_source


def set_classic_plot_style():
    """
    Set the plotting style of the diagnostic figures.

    Returns:
        None
    """
    plt.rcParams.update(
        {
            "font.family": "STIXGeneral",
            "mathtext.fontset": "stix",
            "font.size": 11,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.linewidth": 0.9,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.size": 4,
            "ytick.major.size": 4,
            "xtick.minor.size": 2,
            "ytick.minor.size": 2,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.dpi": 300,
        }
    )


def dist_less_radius(
    flight_utm_x_km,
    flight_utm_y_km,
    seismo_utm_x_km,
    seismo_utm_y_km,
    max_distance_km,
):
    """
    Whether any position of the flight is within the radius of any node.

    Args:
        flight_utm_x_km (numpy.ndarray): UTM x of the flight in km.
        flight_utm_y_km (numpy.ndarray): UTM y of the flight in km.
        seismo_utm_x_km (numpy.ndarray): UTM x of the nodes in km.
        seismo_utm_y_km (numpy.ndarray): UTM y of the nodes in km.
        max_distance_km (float): Radius in km.

    Returns:
        bool: True if the flight comes within the radius of a node.
    """
    for i in range(len(flight_utm_x_km)):
        dx = seismo_utm_x_km - flight_utm_x_km[i]
        dy = seismo_utm_y_km - flight_utm_y_km[i]
        dist_km = np.sqrt(dx * dx + dy * dy)

        if np.any(dist_km <= max_distance_km):
            return True

    return False


_STATION_CACHE = None


def load_station_geometry():
    # read the node file once per worker process, not once per flight
    """
    Read the node file and add UTM coordinates (read once per worker process).

    Returns:
        pandas.DataFrame: One row per node with utm_x_km and utm_y_km.
    """
    global _STATION_CACHE
    if _STATION_CACHE is not None:
        return _STATION_CACHE.copy()
    seismo_data = pd.read_csv(NODE_FILE, sep="|")

    required = ["Latitude", "Longitude", "Station", "StartTime", "EndTime"]
    missing = [col for col in required if col not in seismo_data.columns]

    if missing:
        raise ValueError(f"Missing required columns in {NODE_FILE}: {missing}")

    seismo_utm = [
        utm_proj(lon, lat)
        for lat, lon in zip(seismo_data["Latitude"], seismo_data["Longitude"])
    ]

    seismo_utm_x_m, seismo_utm_y_m = zip(*seismo_utm)

    seismo_data["utm_x_km"] = np.asarray(seismo_utm_x_m, dtype=float) / 1000.0
    seismo_data["utm_y_km"] = np.asarray(seismo_utm_y_m, dtype=float) / 1000.0

    seismo_data["Station"] = seismo_data["Station"].astype(str).str.strip()

    _STATION_CACHE = seismo_data
    return seismo_data.copy()


def safe_get_equip(date, flight_num):
    """
    Aircraft type of a flight, or unknown.

    Args:
        date (str): Date, YYYYMMDD.
        flight_num (str): Flight id.

    Returns:
        str: Equipment code, or unknown.
    """
    try:
        equip = get_equip(date, flight_num)
    except Exception:
        equip = "unknown"

    equip = str(equip).strip()

    if equip == "" or equip.lower() in {"nan", "none", "null"}:
        return "unknown"

    return equip


def read_write_closest(args):
    """
    Compute the crossings of one flight with every node within the radius.

    Args:
        args (tuple): Flight file, file name, task number, number of tasks,
            radius in km, and whether to use the corrected database.

    Returns:
        tuple: Crossing lines (12 columns each) and the crossings whose speed came
        from the corrected positions.
    """
    (
        flight_file,
        filename,
        tracer,
        total_files,
        max_horizontal_distance_km,
        use_corrected,
    ) = args

    if tracer % 100 == 0:
        print(f"{100.0 * tracer / total_files:.2f}%")

    hold_lines = []
    path_speed_rows = []

    seismo_data = load_station_geometry()

    seismo_utm_x_km = seismo_data["utm_x_km"].to_numpy(dtype=float)
    seismo_utm_y_km = seismo_data["utm_y_km"].to_numpy(dtype=float)

    stations = seismo_data["Station"].to_numpy()
    start_times = seismo_data["StartTime"].to_numpy()
    end_times = seismo_data["EndTime"].to_numpy()

    try:
        if use_corrected:
            loaded = load_corrected_flight_file(flight_file, filename)
        else:
            loaded = load_original_flight_file(flight_file, filename)
    except Exception as exc:
        print(f"failed loading flight file {flight_file}: {exc}")
        return hold_lines, path_speed_rows

    if loaded is None:
        return hold_lines, path_speed_rows

    (
        flight_utm_x_km,
        flight_utm_y_km,
        flight_path,
        timestamp,
        alt,
        speed,
        head,
        flight_num,
        date,
        speed_source,
    ) = loaded

    equip = safe_get_equip(date, flight_num)

    if not dist_less_radius(
        flight_utm_x_km=np.asarray(flight_utm_x_km, dtype=float),
        flight_utm_y_km=np.asarray(flight_utm_y_km, dtype=float),
        seismo_utm_x_km=seismo_utm_x_km,
        seismo_utm_y_km=seismo_utm_y_km,
        max_distance_km=max_horizontal_distance_km,
    ):
        return hold_lines, path_speed_rows

    # only nodes that can be within the radius of the flight: a node further than the
    # radius outside the box around the flight cannot be that close to any point of
    # the path, so skipping it changes nothing in the result, it only saves time
    fx = np.asarray(flight_utm_x_km, dtype=float)
    fy = np.asarray(flight_utm_y_km, dtype=float)
    near = np.flatnonzero(
        (seismo_utm_x_km >= np.nanmin(fx) - max_horizontal_distance_km)
        & (seismo_utm_x_km <= np.nanmax(fx) + max_horizontal_distance_km)
        & (seismo_utm_y_km >= np.nanmin(fy) - max_horizontal_distance_km)
        & (seismo_utm_y_km <= np.nanmax(fy) + max_horizontal_distance_km)
    )

    for s in near:
        seismometer = (
            seismo_utm_x_km[s],
            seismo_utm_y_km[s],
        )

        station = stations[s]

        try:
            closest_p, horizontal_distance_km, index = find_closest_point(
                flight_path,
                seismometer,
            )
        except Exception as exc:
            print(f"find_closest_point failed for flight={flight_num}, station={station}: {exc}")
            continue

        if index is None:
            continue

        if horizontal_distance_km > max_horizontal_distance_km:
            continue

        closest_x_km, closest_y_km = closest_p

        try:
            closest_time = closest_time_calc(
                closest_p,
                flight_path,
                timestamp,
                index,
            )
        except Exception as exc:
            print(f"closest_time_calc failed for flight={flight_num}, station={station}: {exc}")
            continue

        if time_check(closest_time, start_times, end_times, s):
            continue

        try:
            alt_avg_m, speed_avg_mps, head_avg = avg_return(
                alt,
                speed,
                head,
                index,
            )
        except Exception as exc:
            print(f"avg_return failed for flight={flight_num}, station={station}: {exc}")
            continue

        horizontal_distance_m = float(horizontal_distance_km) * 1000.0
        closest_x_m = float(closest_x_km) * 1000.0
        closest_y_m = float(closest_y_km) * 1000.0

        slant_distance_m = np.sqrt(horizontal_distance_m**2 + float(alt_avg_m)**2)
        slant_distance_km = slant_distance_m / 1000.0

        line = (
            f"{date},"
            f"{flight_num},"
            f"{closest_x_m},"
            f"{closest_y_m},"
            f"{horizontal_distance_m},"
            f"{closest_time},"
            f"{alt_avg_m},"
            f"{speed_avg_mps},"
            f"{head_avg},"
            f"{station},"
            f"{equip},"
            f"{slant_distance_km}\n"
        )

        hold_lines.append(line)

        # keep a note when the speed came from the corrected positions (reported speed 0, or replaced)
        if speed_source is not None and {speed_source[index], speed_source[index + 1]} & {"path", "vcalc"}:
            path_speed_rows.append((date, flight_num, station, closest_time, speed_avg_mps))

    return hold_lines, path_speed_rows


def make_bulk_crossing_maps(crossing_file, out_fig_dir, max_horizontal_distance_km):
    """
    Make the diagnostic figures of a crossing file.

    Args:
        crossing_file (pathlib.Path): The crossing file.
        out_fig_dir (pathlib.Path): Folder for the figures.
        max_horizontal_distance_km (float): Radius used, for the labels.

    Returns:
        None
    """
    set_classic_plot_style()

    if not crossing_file.exists():
        print(f"cannot make figures because crossing file was not found: {crossing_file}")
        return

    crossings = pd.read_csv(crossing_file, names=CROSSING_COLUMNS, header=None)
    stations = load_station_geometry()

    for col in [
        "X_m",
        "Y_m",
        "HorizontalDistance_m",
        "ClosestTime",
        "Altitude_m",
        "Speed_mps",
        "Heading",
        "SlantDistance_km",
    ]:
        crossings[col] = pd.to_numeric(crossings[col], errors="coerce")

    crossings = crossings.dropna(
        subset=[
            "X_m",
            "Y_m",
            "HorizontalDistance_m",
            "ClosestTime",
            "SlantDistance_km",
            "Speed_mps",
        ]
    ).copy()

    if crossings.empty:
        print("no crossings available for figures")
        return

    crossings["X_km"] = crossings["X_m"] / 1000.0
    crossings["Y_km"] = crossings["Y_m"] / 1000.0
    crossings["HorizontalDistance_km"] = crossings["HorizontalDistance_m"] / 1000.0

    crossings["Date"] = crossings["Date"].astype(str).str.strip()
    crossings["Station"] = crossings["Station"].astype(str).str.strip()

    station_counts = crossings.groupby("Station").size().reset_index(name="crossing_count")
    station_plot = stations.merge(station_counts, on="Station", how="left")
    station_plot["crossing_count"] = station_plot["crossing_count"].fillna(0)

    fig, ax = plt.subplots(figsize=(7.6, 7.0))
    hb = ax.hexbin(crossings["X_km"], crossings["Y_km"], gridsize=70, mincnt=1, bins="log", linewidths=0.0)
    ax.scatter(stations["utm_x_km"], stations["utm_y_km"], s=12, marker="^", facecolors="white",
               edgecolors="black", linewidths=0.35, label="ZE nodes", zorder=5)
    cb = fig.colorbar(hb, ax=ax, pad=0.01, fraction=0.04)
    cb.set_label("Crossing density, log count")
    ax.set_xlabel("UTM Easting (km)")
    ax.set_ylabel("UTM Northing (km)")
    ax.set_title(f"All-aircraft closest-approach density, d0 <= {max_horizontal_distance_km:g} km")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(frameon=False, loc="best")
    ax.minorticks_on()
    fig.tight_layout()
    fig.savefig(out_fig_dir / f"01_crossing_density_map_{radius_label(max_horizontal_distance_km)}km.png",
                dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    ax.hist(crossings["Speed_mps"], bins=np.arange(0, 360, 5), edgecolor="black", linewidth=0.45, alpha=0.85)
    ax.set_xlabel("Speed at crossing (m/s)")
    ax.set_ylabel("Crossing count")
    ax.set_title("Speed at crossing")
    ax.minorticks_on()
    fig.tight_layout()
    fig.savefig(out_fig_dir / f"07_speed_distribution_{radius_label(max_horizontal_distance_km)}km.png",
                dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    daily_counts = crossings.groupby("Date").size().reset_index(name="crossing_count").sort_values("Date")
    daily_counts.to_csv(out_fig_dir / f"daily_crossing_summary_{radius_label(max_horizontal_distance_km)}km.csv",
                        index=False)

    print()
    print("diagnostic figures complete")
    print(f"figures written to: {out_fig_dir}")


def main():
    """
    Compute the crossings of all flights in the date window and write the crossing file.

    Returns:
        None
    """
    args = parse_args()

    max_horizontal_distance_km = float(args.radius_km)
    num_workers = int(args.workers)

    label = radius_label(max_horizontal_distance_km)

    # the original mode gets its own name, so it never overwrites an older crossing file
    tag = "_corrected_db" if args.use_corrected else "_original"
    output_file = Path(f"input/all_station_crossing_db_UTM_{label}km{tag}.txt")
    path_speed_file = Path(f"input/all_station_crossing_db_UTM_{label}km{tag}_path_speed.csv")
    out_fig_dir = Path(f"outputs_aircraft_crossings_{label}km{tag}")
    out_fig_dir.mkdir(parents=True, exist_ok=True)

    if args.use_corrected:
        flight_files, filenames = corrected_flight_list(args.corrected_db, args.start_date, args.end_date)
    else:
        flight_files, filenames = original_flight_list(RAW_DIR, args.start_date, args.end_date)
    print(f"mode: {'corrected database, smoothed speed' if args.use_corrected else 'original FR24 files, reported speed'}")

    total_files = len(flight_files)

    print(f"flights from: {args.corrected_db if args.use_corrected else RAW_DIR}")
    print(f"dates: {args.start_date} to {args.end_date}")
    print(f"flight files: {total_files}")
    print(f"horizontal search radius: {max_horizontal_distance_km:g} km")
    print(f"workers: {num_workers}")
    print(f"output file: {output_file}")

    if total_files == 0:
        print("no flight files found for these dates" + ("; run correct_flight_paths.py first" if args.use_corrected else ""))
        return

    tasks = [
        (
            flight_files[i],
            filenames[i],
            i,
            total_files,
            max_horizontal_distance_km,
            args.use_corrected,
        )
        for i in range(total_files)
    ]

    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        output_list = list(executor.map(read_write_closest, tasks))

    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w") as output:
        for flight_lines, _ in output_list:
            for line in flight_lines:
                output.write(line)

    path_rows = [row for _, rows in output_list for row in rows]
    pd.DataFrame(path_rows, columns=["Date", "FlightNum", "Station", "ClosestTime", "Speed_mps"]).to_csv(
        path_speed_file, index=False
    )

    n_rows = sum(len(flight_lines) for flight_lines, _ in output_list)

    print()
    print("crossing generation complete")
    print(f"rows written: {n_rows}")
    print(f"crossings with speed from the corrected positions: {len(path_rows)}")
    print(f"wrote: {output_file}")
    print(f"wrote: {path_speed_file}")
    print()
    print("columns:")
    print(
        "Date,FlightNum,X_m,Y_m,HorizontalDistance_m,ClosestTime,"
        "Altitude_m,Speed_mps,Heading,Station,Equipment,SlantDistance_km"
    )

    make_bulk_crossing_maps(
        crossing_file=output_file,
        out_fig_dir=out_fig_dir,
        max_horizontal_distance_km=max_horizontal_distance_km,
    )


if __name__ == "__main__":
    main()
