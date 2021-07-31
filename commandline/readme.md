# Command line version

## Arguments

```bash
swEOS -D1 -VT -R0/2/1000 -P316 -X0.032
```

### Required arguments

* -D[0|1|2|3]: calculation model (**D**imension) selection, 0~3 represent scatter point, 1D, 2D, 3D calculation, respectively

* -V[PTX|PHX|T|P|X|H|PT|PX|TX|PH|HX]: **V**ariables selection. **Note** that different `-D` parameter correspond to different options for `-V`.

    1) -D0: [PTX|PHX]
    2) -D1: [T|P|X|H]
    3) -D2: [PT|PX|TX|PH|HX]
    4) -D3: [PTX|PHX]

* -P: fixed pressure value

* -T: fixed temperature value

* -X: fixed salinity value

* -H: fixed enthalpy value

* -R: only for 1~3 dimension. Corresponding to `-V` parameter, e.g. `-VPT` and `-R10/1/600/0/1/1000` means pressure in range of [10, 600] bar and interval is 1 bar, temperature in range of [0, 1000] deg. C and interval is 1 deg. C.

* -G: input file for multiple points calculation

* -O: output file for 1D, 2D, and 3D calculation

### Optional arguments

* -t: number of threads




## Example
### Single point calculation
```
swEOS -D0 -VPXT -P316 -T100 -X0.032
swEOS -D0 -VPXH -H438 -P316 -X0.032
```
### Multi-points calculation
```
swEOS -D0 -VPHX -G../test/PHX.txt -OPHX_0D.csv
swEOS -D0 -VPHX -G../test/PHX.txt
swEOS -D0 -VPTX -G../test/PTX.txt -OPTX_0D.csv
swEOS -D0 -VPTX -G../test/PTX.txt
```

### One-dimensional calculation

```
swEOS -D1 -VH -X0.032 -P399 -R43/1/100 -OH_1D.csv
swEOS -D1 -VT -X0.032 -P399 -R0/1/100 -OT_1D.csv
swEOS -D1 -VP -X0.032 -T100 -R5/1/100 -OP_1D.csv
swEOS -D1 -VX -T100 -P399 -R0/0.001/1 -OX_1D.csv
```

### Two-dimensional calculation

```
swEOS -D2 -VPT -R1/0.1/100/0/1/500 -X0.032 -Orest/PT_2D.vt 
swEOS -D2 -VPX -R100/0.1/800/0/0.01/1 -T100  -OPX_2D.vtk -n
swEOS -D2 -VTX -R0/1/800/0/0.01/1 -P100  -OTX_2D.vtk -n
swEOS -D2 -VPH -R100/1/800/100/1/700 -X0.032  -OPH_2D.vtk -n
swEOS -D2 -VXH -R0/0.001/1/100/1/700 -P200  -OXH_2D.vtk -n
```

### Three-dimensional calculation

```
swEOS -D3 -VPTX -R1/10/500/0/10/600/0/0.01/1 -t8 -n
swEOS -D3 -VPHX -R1/10/500/100/10/600/0/0.01/1 -t8 -n
```