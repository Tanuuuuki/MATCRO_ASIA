##### SET for MATCRO #####
touch ./RUN/MATCRO.f90
cd ./RUN; make ; wait; cd ..
echo "make OK"

country=$1
if [ ${#country} -ne 3 ]; then
    echo "wrong argument"
    echo $1
    exit
fi

opt_program=./PROGRAM/double_optim_LCB.py
orifile=./SET/SET_CIRES_ori.txt

for region in $2 ; do    # Set regionecture for optimization
    
    list_dir=./$country/$region
    mkdir -p $list_dir
    list_file=$list_dir/${KOR}_${region}.list
    ls ./CLM/tmp/$country/* > $list_file

    ## decide year_st and year_ed ##
    stat_file=./STATS/${country}_yield.csv
    year_st=`grep $region $stat_file | cut -d',' -f3 | head -n 1`
    year_ed=`grep $region $stat_file | cut -d',' -f3 | tail -n 1`
    echo $year_st $year_ed
    if [ $year_st -lt 1896 ]; then year_st=1896; fi
    if [ $year_ed -gt 1996 ]; then year_ed=1996; fi
    year_list=(`seq $year_st 10 $year_ed | tr '\n' ' '`)
    echo $year_list

    ## year reversing to start from recent years ##
    for ((i=${#year_list[@]}-1; i>=0; i--)); do
        reversed_year+=("${year_list[i]}")
    done

    ##### DISTRICT CSV row #####
    csv_row=`grep $region ~/CRU/DISTRICT_AREA/code_name_0.5.csv | cut -d',' -f1`

    #echo $reversed_year
    #exit 1
    for year in $reversed_year
    do
        echo $year

        opt_term=20                                          # Set optimization period
        styear=$year                                         # start year
        edyear=$(($styear+$opt_term-1))                      # end year
        year_grep=`seq $styear 1 $edyear | tr '\n' ' ' | sed -E "s/.$/\n/" | sed -e "s/ /|/g"`
        #echo "($year_grep)"

        ##### set lonlat #####
        longitude=`grep $year ~/MATCRO/Baye2/INPUT/Coordinate/Sap_lonlat.csv | cut -d',' -f2`
        latitude=`grep $year ~/MATCRO/Baye2/INPUT/Coordinate/Sap_lonlat.csv | cut -d',' -f3`        

        season_flag=0
        for season in 1st 2nd 3rd; do
            check_season=`grep $season $list_file | grep -E "($year_grep)" | wc -l`
            #grep $region $clm_file | grep $season | grep | grep -E "($year_grep)"

            if [ ${check_season} -gt 0 ]; then
                echo $season
                (( season_flag++ ))

                check_ir=`grep $season $list_file | grep -E "($year_grep)" | grep ir | wc -l`
                check_rf=`grep $season $list_file | grep -E "($year_grep)" | grep rf | wc -l`
                echo $check_ir $check_rf

                inifile_ir=./SET/$country/$region/SET_${region}_ir${season}_INI.txt
                inifile_rf=./SET/$country/$region/SET_${region}_rf${season}_INI.txt
                optfile_ir=./SET/$country/$region/SET_${region}_ir${season}_opt.txt
                optfile_rf=./SET/$country/$region/SET_${region}_rf${season}_opt.txt

                ## file setting ##
                ##### make directory #####
                mkdir -p ./SET/$country/$region
                mkdir -p ./SMPL/$country/$region
                mkdir -p ./YLD/$country/$region
                mkdir -p ./PRM/$country/$region

                ##### create regionectural setting file for initial and optim #####
                sed "s/Country/$country/g" $orifile | sed "s/Region/$region/g" |\
                sed "/STYR/s/10000$/$styear/" | sed "/ENYR/s/20000$/$edyear/" |\
                sed "/PREFLON/s/lon$/$longitude/g" | sed "/PREFLAT/s/lat$/$latitude/g" |\
                sed "s/season/$season/g" | sed "/IRRF_NAME/s/irrf/ir/" | sed "/IRR /s/irrf/1/" |\
                sed "s/irrf_/ir/g" | sed 's/ori/INI/g' | sed "/PLT_ROW/s/csv_row/$csv_row/" \
                > $inifile_ir

                sed "s/Country/$country/g" $orifile | sed "s/Region/$region/g" |\
                sed "/STYR/s/10000$/$styear/" | sed "/ENYR/s/20000$/$edyear/" |\
                sed "/PREFLON/s/lon$/$longitude/g" | sed "/PREFLAT/s/lat$/$latitude/g" |\
                sed "s/season/$season/g" | sed "/IRRF_NAME/s/irrf/ir/" | sed "/IRR /s/irrf/1/"|\
                sed "s/irrf_/ir/g" | sed 's/ori/opt/g' | sed "/PLT_ROW/s/csv_row/$csv_row/" \
                > $optfile_ir

                sed "s/Country/$country/g" $orifile | sed "s/Region/$region/g" |\
                sed "/STYR/s/10000$/$styear/" | sed "/ENYR/s/20000$/$edyear/" |\
                sed "/PREFLON/s/lon$/$longitude/g" | sed "/PREFLAT/s/lat$/$latitude/g" |\
                sed "s/season/$season/g" | sed "/IRRF_NAME/s/irrf/rf/" | sed "/IRR /s/irrf/0/" |\
                sed "s/irrf_/rf/g" | sed 's/ori/INI/g' | sed "/PLT_ROW/s/csv_row/$csv_row/" \
                > $inifile_rf

                sed "s/Country/$country/g" $orifile | sed "s/Region/$region/g" |\
                sed "/STYR/s/10000$/$styear/" | sed "/ENYR/s/20000$/$edyear/" |\
                sed "/PREFLON/s/lon$/$longitude/g" | sed "/PREFLAT/s/lat$/$latitude/g" |\
                sed "s/season/$season/g" | sed "/IRRF_NAME/s/irrf/rf/" | sed "/IRR /s/irrf/0/" |\
                sed "s/irrf_/rf/g" | sed 's/ori/opt/g' | sed "/PLT_ROW/s/csv_row/$csv_row/" \
                > $optfile_rf

            fi
        done
        echo season_flag is $season_flag


        ##### READ AND WRITE PARAMETER FILE #####
        ### z = zero , f = first , s = second , t = third
        path=./OUTPUT/$country/$region       # regionecture path
        zodir=${path}/BFGS_1-0                 # output path
        zoydir=${zodir}/OPT_YLD_               # optimized yield path
        zopdir=${zodir}/OPT_PRM_               # optimized prm path

        fodir=${path}/BFGS_1-1                 # output path
        foydir=${fodir}/OPT_YLD_               # optimized yield path
        fopdir=${fodir}/OPT_PRM_               # optimized prm path

        sodir=${path}/BFGS_1-2                 # output path
        soydir=${sodir}/OPT_YLD_               # optimized yield path
        sopdir=${sodir}/OPT_PRM_               # optimized prm path

        todir=${path}/BFGS_1                   # output path
        toydir=${todir}/OPT_YLD_               # optimized yield path
        topdir=${todir}/OPT_PRM_               # optimized prm path

        ### mkdir ###
        mkdir -p $fodir $sodir $todir #$zodir
        
        echo optimization start
        ##### first optimization for SLNYMN HI mGDD #####
        #python3 $opt_program \
        #    $country $region $styear $edyear $zodir $topdir $zoydir $zopdir "zero" $season_flag \
        #    | tee ${zodir}/${styear}-${edyear}.out
        #wait

        #echo "wait for first optim"
        #sleep 5

        python3 $opt_program \
            $country $region $styear $edyear $fodir $zopdir $foydir $fopdir "first" $season_flag \
            | tee ${fodir}/${styear}-${edyear}.out
        wait

        echo "wait for second optim"
        sleep 5

        ##### second optimization for TCmin #####
        python3 $opt_program \
            $country $region $styear $edyear $sodir $fopdir $soydir $sopdir "second" $season_flag \
            | tee ${sodir}/${styear}-${edyear}.out
        wait

        ##### final optimization for phenology #####
        python3 $opt_program \
            $country $region $styear $edyear $todir $sopdir $toydir $topdir "third" $season_flag \
            | tee ${todir}/${styear}-${edyear}.out
        wait

        sed -i "/PREFLON/s/[0-9].*[0-9]$/lon/g" $inifile_ir $inifile_rf $optfile_ir $optfile_rf
        sed -i "/PREFLAT/s/[0-9].*[0-9]$/lat/g" $inifile_ir $inifile_rf $optfile_ir $optfile_rf
        sed -i "/STYR/s/[12][0-9][0-9]6$/10000/g" $inifile_ir $inifile_rf $optfile_ir $optfile_rf
        sed -i "/ENYR/s/[12][0-9][0-9]5$/20000/g" $inifile_ir $inifile_rf $optfile_ir $optfile_rf

        grep '\(ST\|EN\)YR' $inifile_ir
        grep 'PREF\(LON\|LAT\)' $optfile_rf

        exit 1        

    done
done

