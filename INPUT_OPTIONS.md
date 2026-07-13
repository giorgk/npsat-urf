# NPSAT URF input options

The C++ executable reads options from `npsat_urf.opt` in the working directory where the binary is run. Each option is read as:

```text
option_name value
```

The parser reads the first token as the option name and then reads the next value with C++ stream extraction, so comments after the value are ignored naturally.

## Input and output files

| Option | Aliases | Default | How it is used |
| --- | --- | --- | --- |
| `prefix` | none | none | Prefix for the input streamline file. The full input filename is `prefix` + zero-padded process id + `.` + `suffix`. |
| `paddingZeros` | none | uninitialized | Number of digits used by the zero-padded process id in the input filename. |
| `suffix` | none | none | Input filename suffix after the dot. |
| `file_type` | none | `npsat_ascii` | Input streamline format. Supported values are `npsat_ascii` and `npsat_bin`; `modpath` is recognized but not implemented. |
| `output_prefix` | none | none | Prefix for the URF output file. The output filename is `output_prefix_<process_id>.dat`. |
| `discard_prefix` | `discarded_prefix`, `discard_output_prefix` | `discarded_streamlines` | Prefix for discarded-streamline diagnostics. The file is `discard_prefix_<process_id>.dat`. |
| `simplified_prefix` | `simplify_prefix` | `simplified_streamline` | Prefix for simplified streamline outputs. The CSV file is `simplified_prefix_<process_id>.dat`; the VTK file is `simplified_prefix_<process_id>.vtk`. |

## Streamline simplification and VTK

| Option | Aliases | Default | How it is used |
| --- | --- | --- | --- |
| `simplify_streamline` | `simplifyStreamline` | `0` | When nonzero, writes simplified streamline points to `simplified_prefix_<process_id>.dat`. |
| `simplify_tolerance` | `simplification_tolerance` | `0.0` | Douglas-Peucker distance tolerance used to simplify each streamline. A value of `0.0` keeps all geometrically non-collinear points. |
| `write_simplified_vtk` | `simplified_vtk`, `writeSimplifiedVtk` | `0` | When nonzero, also writes `simplified_prefix_<process_id>.vtk` for ParaView. This automatically enables `simplify_streamline`. |

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

