#!/bin/bash

#sort files
#for a in CRI* ; do set=$(dfits $a | fitsort -d INS.WLEN.ID | cut -f2 ); mkdir -p $set; mv $a $set; done

mksetting () {
    echo $BASE
    s=$1
    cd $s
    rm -rf *_out *.sof
    cp -r $BASE/sofs/* .

    cp /home/tom/pipes/cr2rep.git/catalogs/selections/$s.dat une_sel.dat

    echo "../cr2res_cal_detlin_coeffs.fits CAL_DETLIN_COEFFS" | tee -a 02_cr2res_util_calib.sof | tee -a 07_cr2res_util_calib.sof >> 13_cr2res_util_calib.sof

    dfits CR*fits | fitsort OBJECT | grep DARK > 01_cr2res_cal_dark.sof
    dfits CR*fits | fitsort OBJECT | grep FLAT >> 02_cr2res_util_calib.sof
    dfits CR*fits | fitsort OBJECT | grep "WAVE,UNE"  | sed s/,UNE/_UNE/ >> 07_cr2res_util_calib.sof
    dfits CR*fits | fitsort OBJECT | grep "WAVE,FPET" | sed s/,FPET/_FPET/ >> 13_cr2res_util_calib.sof

    tail -1 13_cr2res_util_calib.sof > 04_cr2res_util_slit_curv.sof
    echo "03_cr2res_util_trace_out/cr2res_util_calib_calibrated_collapsed_tw.fits  CAL_FLAT_TW" >> 04_cr2res_util_slit_curv.sof

    $BASE/getDarkname_noNDIT.py `tail -1 02_cr2res_util_calib.sof| awk '{print $1}'` | awk '{print  "01_cr2res_cal_dark_out/" $0}' | sponge >> 02_cr2res_util_calib.sof
    $BASE/getDarkname_noNDIT.py `tail -1 07_cr2res_util_calib.sof| awk '{print $1}'` | awk '{print  "01_cr2res_cal_dark_out/" $0}' | sponge >> 07_cr2res_util_calib.sof
    $BASE/getDarkname_noNDIT.py `tail -1 13_cr2res_util_calib.sof| awk '{print $1}'` | awk '{print  "01_cr2res_cal_dark_out/" $0}' | sponge >> 13_cr2res_util_calib.sof

    echo "06_cr2res_util_normflat_out/cr2res_util_normflat_Open_master_flat.fits CAL_FLAT_MASTER" >> 07_cr2res_util_calib.sof
    echo "06_cr2res_util_normflat_out/cr2res_util_normflat_Open_master_flat.fits CAL_FLAT_MASTER" >> 13_cr2res_util_calib.sof

    if [[ $s == M* ]] || [[ $s == L* ]] ; then
     replace 04_cr2res_util_slit_curv_out/cr2res_util_calib_calibrated_collapsed_tw_tw.fits 03_cr2res_util_trace_out/cr2res_util_calib_calibrated_collapsed_tw.fits -- *sof
    fi
    if [[ $s == Y* ]] || [[ $s == J* ]] ; then
     replace redman sarmiento -- *sof
    fi

}
export -f mksetting


BASE=`readlink -f /home/tom/pCOMM`
cp -s $BASE/210313_detlin/cr2res_cal_detlin_coeffs.fits .
cp /home/tom/pipes/cr2rep.git/catalogs/lines_u* .

parallel "BASE=$BASE mksetting {}" ::: [YJHKLM][1-4][0-9][0-9][0-9]


exit 0;
