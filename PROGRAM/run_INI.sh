#!bin/nash
cd ./RUN/INI/; make; cd ../../

BFGS=$1
if [ $# = 0 ]; then
    echo you have to set BFGS file
    exit 1
fi

country=KOR 		# south Korea
code_name_path=~/CRU/DISTRICT_AREA/code_name_0.5.csv
echo $country
echo $code_name_path
result=`grep $country $code_name_path | cut -d"," -f3 | tr '\n' ' ' | xargs echo`
IFS=' ' read -r -a district <<< $result
dir_ann=(YLD PRM LAImx WSH WSO WAR WST WLF WRT WDL)
dir_day=(LAI SLN GPP TSP RSP)

#for name in ${district[@]}; do echo $name; done
#exit 1

##### make initial SET file #####
orifile=./SET/SET_CIRES_ori.txt
modfile=./SET/SET_CIRES_initial.txt
sed "s/ori/initial/g" $orifile > $modfile

  styear=1896
  edyear=1915

#for region in $@    # Set regionecture for initialization
for region in ${district[@]}; do # "Chungcheongnam-do" "Gyeonggi-do" "Gyeongsangnam-do"; do
  echo $region

  out_path=./ANALYSIS/$BFGS/$country/$region
  for DIR in ${dir_ann[@]} ; do
    out_dir=$out_path/$DIR
    tmp_dir=$DIR/$country/$region
    mkdir -p $out_dir $tmp_dir
  done

  for DIR in ${dir_day[@]} ; do
    out_dir=$out_path/OPT_daily/$DIR
    tmp_dir=$DIR/$country/$region
    mkdir -p $out_dir $tmp_dir
  done


  ## regionecture setting ##
  tmpfile=./SET/$country/$region/SET_${region}_tmp.txt
  
  ##### create regionectural setting file #####
  mkdir -p ./SET/$country/$region

  sed "s/Region/$region/g" $modfile > $tmpfile
  sed -i "s/Country/$country/g" $tmpfile

  ##### set STYR and ENYR in setfile #####
  sed -i "s/10000$/$styear/g" $tmpfile
  sed -i "s/20000$/$edyear/g" $tmpfile

  ##### irrigation or rainfed #####
  for irrf in ir rf; do
    echo $irrf

    if [ $irrf = ir ]; then
      sed -i "/IRR /s/irrf/1/g" $tmpfile
    elif [ $irrf = rf ]; then
      sed -i "/IRR /s/irrf/0/g" $tmpfile
    else
      echo "something wrong happend"
      exit 1
    fi
    sed -i "s/irrf/$irrf/g" $tmpfile

    grep IRR $tmpfile
      
    check_irrf=./CLM/tmp/$country/${styear}_${irrf}1st_tmp_${region}.txt
    echo $check_irrf

    if [ -e $check_irrf ]; then

      ##### RUN for 1st, 2nd, and 3rd #####
      for season in 1st 2nd 3rd; do
        echo $season
        sed -i "s/season/$season/g" $tmpfile
        grep SEASON $tmpfile

        check_season=./CLM/tmp/$country/${styear}_${irrf}${season}_tmp_${region}.txt
        if [ -e $check_season ]; then
            ./RUN/INI/MATCRO_INI $tmpfile
      
          for DIR in ${dir_ann[@]} ; do
            ini_file=./$DIR/$country/$region/${DIR}_${region}_${irrf}_${season}_initial.txt
            out_file=$out_path/$DIR/${DIR}_${region}_${irrf}_${season}_initial.txt
            if [ ! -e $ini_file ]; then
              echo "error happend"
              echo $ini_file $out_file
              exit 1
            fi

            if [ $DIR == "YLD" ]; then echo $out_file; fi
            mv $ini_file $out_file
          done

          for DIR in ${dir_day[@]} ; do
            out_dir=$out_path/OPT_daily/$DIR

            for tyear in $(seq $styear $edyear); do
              ifile=./$DIR/$country/$region/${DIR}_${region}_${tyear}.csv
              if [ ! -e $ifile ]; then
                echo "error happend"
                echo $styear $edyear $tyear
                echo $ifile
                exit 1
              fi

              ofile=$out_dir/${DIR}_${region}_${tyear}_initial.csv
              mv $ifile $ofile
            done
            
          done  
        fi #if season clm file exist

        sed -i "s/$season/season/g" $tmpfile
        #grep SEASON $tmpfile

      done #season
    fi

    sed -i "s/$irrf/irrf/g" $tmpfile

    if [ $irrf = ir ]; then
      sed -i "/IRR /s/1/irrf/g" $tmpfile
    elif [ $irrf = rf ]; then
      sed -i "/IRR /s/0/irrf/g" $tmpfile
    fi

  done #irrf

done #region

##### we have to add command to calculate total yld
## Rscript @@