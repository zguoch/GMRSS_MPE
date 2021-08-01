gmt begin wave_t pdf 
    gmt basemap -R0/500/-7/7 -JX12c/6c -Ba 
    awk '{print $3, $4}' test_h.txt  | gmt plot -W0.2p
gmt end show