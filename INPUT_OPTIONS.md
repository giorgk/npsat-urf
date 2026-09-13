# NPSAT URF input options

The C++ executable reads options from the file passed as its sole command-line argument. The filename is not fixed; `npsat_urf.opt` and `option_file.dat` are both valid names. Relative input and output paths are resolved from the working directory where the binary is run. Each option is read as:

```text
option_name value
```

The parser reads the first token as the option name and then reads the next value with C++ stream extraction, so comments after the value are ignored naturally.

## Input and output files

| Option | Aliases | Default | How it is used |
| --- | --- | --- | --- |
| `prefix` | none | none | Prefix before the built-in `rank_<rank>_iter_<iter>` filename segment. |
| `paddingZeros` | none | uninitialized | Number of digits used by the zero-padded rank id in the input filename. |
| `iter_paddingZeros` | `iterPaddingZeros` | `4` | Number of digits used by the zero-padded iteration id in the input filename. |
| `suffix` | none | none | Input filename suffix after the dot. |
| `file_type` | none | `npsat_ascii` | Input streamline format. Supported values are `npsat_ascii` and `npsat_bin`; `modpath` is recognized but not implemented. |
| `output_prefix` | none | none | Prefix for the URF output file. The output filename is `output_prefix_rank_<mpi_rank>.dat`. |
| `discard_prefix` | `discarded_prefix`, `discard_output_prefix` | `discarded_streamlines` | Prefix for discarded-streamline diagnostics. The file is `discard_prefix_rank_<mpi_rank>.dat`. |
| `simplified_prefix` | `simplify_prefix` | `simplified_streamline` | Prefix for simplified streamline outputs. The CSV file is `simplified_prefix_rank_<mpi_rank>.dat`; the VTK file is `simplified_prefix_rank_<mpi_rank>.vtk`. |

## Streamline simplification and VTK

| Option | Aliases | Default | How it is used |
| --- | --- | --- | --- |
| `simplify_streamline` | `simplifyStreamline` | `0` | When nonzero, writes simplified streamline points to `simplified_prefix_rank_<mpi_rank>.dat`. |
| `simplify_tolerance` | `simplification_tolerance` | `0.0` | Douglas-Peucker distance tolerance used to simplify each streamline. A value of `0.0` keeps all geometrically non-collinear points. |
| `write_simplified_vtk` | `simplified_vtk`, `writeSimplifiedVtk` | `0` | When nonzero, also writes `simplified_prefix_rank_<mpi_rank>.vtk` for ParaView. This automatically enables `simplify_streamline`. |

The simplified CSV columns are `Eid,Sid,x,y,z,v,a`. The VTK output is legacy ASCII `POLYDATA`: all simplified points are written once, and each contiguous `Eid/Sid` group is written as one polyline in the `LINES` section. Point data includes velocity, age, `Eid`, and `Sid`.

## URF calculation controls

| Option | Default | How it is used |
| --- | --- | --- |
| `alpha` | `0.32` | Longitudinal dispersivity coefficient used in URF calculations. |
| `beta` | `0.83` | Exponent used with `alpha` in URF calculations. |
| `Dm` | `1.1578e-4` | Molecular diffusion coefficient used when diffusion/decay calculations are enabled. |
| `minElemSize` | `0.01` | Minimum streamline segment length used when building URF calculation segments. Short segments are carried forward until the accumulated length is large enough. |
| `maxElemSize` | `20` | Maximum streamline segment length. Longer segments are split into smaller pieces. |
| `TimeStep` | `365.0` | Time step used in URF calculations. |
| `maxTotalTime` | `1000` | Maximum total time considered by the URF calculation. |
| `URFtol` | `0.99` | URF cumulative tolerance used by the calculation. |
| `skipAge` | `2` | Age-sampling skip parameter used by the URF algorithm. |
| `halfTime` | `12.32` | Half-time parameter used for decay calculations. |
| `er_to_run` | `1` | End reason to process. Values below zero process all end reasons; otherwise only streamlines with matching end reason get URF values. |

`er_to_run` may contain multiple values on the same line. For example:

```text
er_to_run 7 14 3
```

With this setting, fitting is carried out when the streamline end reason is `7`, `14`, or `3`. If any listed value is negative, all end reasons are processed.

## Command-line mode

The current C++ implementation expects exactly one argument: the options-file path.

```bash
mpirun -n <urf_processes> NPSAT_URF <options_file>
```

For a single process, use `NPSAT_URF <options_file>`. Use `NPSAT_URF -v` to print the version.

The previously documented `NPSAT_URF <n_ranks> <n_iters> <array_id>` syntax is obsolete and is rejected by the current code. `array_id` represented a flattened index identifying one input rank/iteration pair. You do not supply it now: the program generates all work indices internally and distributes them across its MPI processes.

### Input layout and execution options

| Option | Default | How it is used |
| --- | --- | --- |
| `nproc` | `0` (must be set) | Number of ranks used by the original `npsat_trace` run, determining the input rank ids `0` through `nproc - 1`. This is independent of the number of URF MPI processes. |
| `niter` | `0` (must be set) | Number of input iterations, determining iteration ids `0` through `niter - 1`. This is a count, not the last iteration id. |
| `output_files_per_part` | `0` | Number of input files processed by each URF MPI process before starting another output part. `0` keeps all that process's results in one part; must be nonnegative. |
| `progress_percent` | `5` | Progress reporting interval in percent for completed files overall and bytes processed within each input file; must be between `1` and `100`. |

Each MPI process prints per-file progress after finishing streamlines, including its rank, work id, input filename, and processed streamline count. The per-file percentage uses the input-file position in bytes; it is not a time estimate. A single expensive streamline can still cause a pause between updates. Overall `Progress` lines count completed input files. If processing crosses several percentage thresholds at once, only the current percentage is printed, rather than replaying the skipped milestones.

The program expects all `nproc * niter` rank/iteration files. For each internal `work_id` from `0` through `nproc * niter - 1`, it computes:

```text
input_rank = work_id % nproc
input_iter = work_id / nproc     (integer division)
```

It reads the corresponding filename:

```text
prefix + "rank_" + zero_pad(input_rank, paddingZeros) + "_iter_" + zero_pad(input_iter, iter_paddingZeros) + "." + suffix
```

URF MPI process `r` handles work ids `r`, `r + urf_processes`, `r + 2 * urf_processes`, and so on, while they remain below `nproc * niter`.

### Example: 16 trace processors, 5 iterations, 6 URF processors

Assume the trace files use input rank ids `0`–`15` and iteration ids `0`–`4`, with names such as `mcm_vi_streamlines_ordered_rank_0000_iter_0000.bin` through `mcm_vi_streamlines_ordered_rank_0015_iter_0004.bin` (80 files total).

Put the following in `option_file.dat`, along with any desired calculation options:

```text
nproc 16
niter 5
prefix mcm_vi_streamlines_ordered_
paddingZeros 4
iter_paddingZeros 4
suffix bin
file_type npsat_bin
output_prefix urf
discard_prefix discarded_streamlines
output_files_per_part 0
progress_percent 5
```

Adjust `prefix`, `suffix`, and `file_type` to match your actual input files. Then run:

```bash
mpirun -n 6 NPSAT_URF option_file.dat
```

Nothing follows `option_file.dat`. The `6` selects the number of URF MPI processes; `nproc 16` and `niter 5` describe the existing trace files. Do not change `nproc` to `6` for this run.

For example, URF MPI process `0` handles work ids `0, 6, 12, 18, ... , 78`. Work id `18` selects input rank `2`, iteration `1`. Processes `0` and `1` each handle 14 files; processes `2`–`5` each handle 13 files.

With `output_files_per_part 0`, the URF outputs are `urf_rank_0.dat` through `urf_rank_5.dat`, each containing results from that URF process's assigned input files. Discard diagnostics are `discarded_streamlines_rank_0.dat` through `discarded_streamlines_rank_5.dat`. Output rank ids refer to the six URF processes, not the original 16 trace ranks.

If `output_files_per_part` is positive, the first part keeps the same filename, and later parts add `_1`, `_2`, etc. before the extension (for example, `urf_rank_0_1.dat`). The same naming rule applies to discard diagnostics and optional simplified CSV/VTK outputs.

## Porosity sweep

| Option | Default | How it is used |
| --- | --- | --- |
| `startPor` | `10` | First porosity multiplier value in the output sweep. The code uses `value / 10.0` as the velocity multiplier. |
| `endPor` | `60` | Last porosity multiplier value in the output sweep. |
| `intervPor` | `10` | Increment between porosity multiplier values. |

## Optional calculation switches

| Option | Default | How it is used |
| --- | --- | --- |
| `bIsGather` | `0` | Parsed and stored for compatibility with existing option files. |
| `calcDecay` | `0` | When nonzero, writes decay fit columns in addition to the base URF columns. |
| `calcDiff` | `0` | When nonzero, writes diffusion fit columns. If `calcDecay` is `0`, this is forced to `0`. |
