#!bin/nash
cd ./RUN/INI/; make; cd ../../

BFGS=$1
if [ $# = 0 ]; then
    echo you have to set BFGS file
    exit 1
fi

country_list=(KOR)
#district=()
dir_ann=(YLD PRM LAImx WSH WSO WAR WST WLF WRT WDL)
dir_day=(LAI SLN GPP TSP RSP)

code_name_path=~/CRU/DISTRICT_AREA/code_name_0.5.csv

orifile=./SET/SET_CIRES_ori.txt
modfile=./SET/SET_CIRES_initial.txt

styear=1896
edyear=2015

for country in ${country_list[@]}; do
  echo $country

  ##### make district list #####
  district_ori=$(grep $country $code_name_path | cut -d"," -f3 | tr -d ' ' | tr '\n' ' ' | xargs echo)
  #IFS=' ' read -r -a district_ori <<< $result
  #echo $result
  echo $district_ori
  #file_exist=`grep $region $check | wc`

  district=()
  for region in $district_ori; do
    file_count=`ls ./CLM/tmp/$country/* | grep $region | wc -l`
    if [ $file_count -gt 0 ]; then
        echo $region
        district+=($region)
        #echo ${district[@]}
        #exit
    fi
  done
  
  echo ${district[@]}
  #exit 1
  #echo $result
  ##### make initial SET file #####
  sed "s/ori/initial/g" $orifile > $modfile

  for region in ${district[@]}; do
    echo $region

    ##### irrigation or rainfed #####
    for irrf in ir rf; do
      echo $irrf

      ##### make output directory ##### ここで初めて設定するのでよい。
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
          
      ##### RUN for 1st, 2nd, and 3rd #####
      for season in 1st 2nd 3rd; do
        echo $season
        sed -i "s/season/$season/g" $tmpfile
        grep SEASON $tmpfile

        #check_season=./CLM/tmp/$country/${styear}_${irrf}${season}_tmp_${region}.txt
        #if [ -e $check_season ]; then
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
            ifile_clm=./CLM/tmp/$country/${tyear}_${irrf}${season}_tmp_${region}.txt
            if [ ! -e $ifile -a ! -e $ifile_clm ]; then
              echo "no output is OK"
              continue
            elif [ ! -e $ifile ]; then
              echo "error happend"
              echo $styear $edyear $tyear
              echo $ifile
              exit 1
            fi

            ofile=$out_dir/${DIR}_${region}_${tyear}_initial.csv
            mv $ifile $ofile
          done
          
        done  
        #fi #if season clm file exist

        sed -i "s/$season/season/g" $tmpfile
        #grep SEASON $tmpfile

      done #season

      sed -i "s/$irrf/irrf/g" $tmpfile

      if [ $irrf = ir ]; then
        sed -i "/IRR /s/1/irrf/g" $tmpfile
      elif [ $irrf = rf ]; then
        sed -i "/IRR /s/0/irrf/g" $tmpfile
      fi
    done #irrf
    
  done #region

  exit
  ##### we have to add command to calculate total yld
  ## Rscript @@

done #country
