# rm test_h.txt 
# ./GMRSSofMPE OceanWave  -R 0/1/0/0/1/0/0/5/10/0/1/500 -M 0.01/3 -U 20 -A 45 -T 1 -t 8 -O test -f txt

# rm test2_h.txt 
# ./GMRSSofMPE OceanWave  -R 0/1/0/0/1/0/0/5/10/0/1/500 -M 0.01/3 -U 20 -A 45 -T 1 -t 8 -O test2 -f txt
gmt gmtset FONT_LABEL=8p,Helvetica,black
gmt begin wave_t pdf 
    gmt basemap -R0/500/-10/10 -JX12c/4c -Ba -Bxa+l"Time (s)" -Bya+l"Ocean Wave Height (m)"
    awk '{print $3, $4}' test_h.txt  | gmt plot -W0.2p,blue
    # awk '{print $3, $4}' test2_h.txt  | gmt plot -W0.2p,red
    gmt basemap -R0/500/-4/4 -JX12c/4c -Ba -Bxa+l"Time (s)" -Bya+l"Velocity (m/s)" -Y-5.5c
    awk '{print $4, $5}' test_U.txt  | gmt plot -W0.2p,blue -l"Ux"
    awk '{print $4, $6}' test_U.txt  | gmt plot -W0.2p,red -l"Uy"
    awk '{print $4, $7}' test_U.txt  | gmt plot -W0.2p,green -l"Uz"
    echo "L 12p" | gmt legend -DjTR+w3c+o0.1c -F+p0.2p+gwhite@10 
gmt end show