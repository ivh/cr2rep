#!/bin/bash
SETT=$1

cd $SETT

# fake GOTO fi
#if false; then

#fi # GOTO
echo "Next ${SETT}: cr2res_cal_dark 01_cr2res_cal_dark.sof"
esorex --log-file=01_cr2res_cal_dark.log --output-dir=01_cr2res_cal_dark_out cr2res_cal_dark 01_cr2res_cal_dark.sof
if [ -e 01_cr2res_cal_dark_out/cr2res_cal_dark_master.fits ] ; then
  newname=$(../../getDarkname.py 01_cr2res_cal_dark_out/cr2res_cal_dark_master.fits | grep MASTER | cut -f1 -d\  )
  mv 01_cr2res_cal_dark_out/cr2res_cal_dark_master.fits "01_cr2res_cal_dark_out/$newname"
fi
if [ -e 01_cr2res_cal_dark_out/cr2res_cal_dark_bpm.fits ] ; then
  newname=$(../../getDarkname.py 01_cr2res_cal_dark_out/cr2res_cal_dark_bpm.fits | grep BPM | cut -f1 -d\  )
  mv 01_cr2res_cal_dark_out/cr2res_cal_dark_bpm.fits "01_cr2res_cal_dark_out/$newname"
fi
find 01_cr2res_cal_dark_out/ -type l | xargs rm -f
mcp -s "01_cr2res_cal_dark_out/cr2res_cal_dark_*_*x*_*.fits" "01_cr2res_cal_dark_out/cr2res_cal_dark_#1_#2_#4.fits"
echo "Next ${SETT}: cr2res_util_calib 02_cr2res_util_calib.sof (calib+comb flats)"
esorex --log-file=02_cr2res_util_calib.log --output-dir=02_cr2res_util_calib_out cr2res_util_calib --collapse="MEAN" 02_cr2res_util_calib.sof
echo "Next ${SETT}: cr2res_util_trace 03_cr2res_util_trace.sof"
esorex --log-file=03_cr2res_util_trace.log --output-dir=03_cr2res_util_trace_out cr2res_util_trace 03_cr2res_util_trace.sof
echo "Next ${SETT}: cr2res_util_slit_curv 04_cr2res_util_slit_curv.sof"
esorex --log-file=04_cr2res_util_slit_curv.log --output-dir=04_cr2res_util_slit_curv_out cr2res_util_slit_curv 04_cr2res_util_slit_curv.sof
echo "Next ${SETT}: cr2res_util_extract 05_cr2res_util_extract.sof (extact flat)"
esorex --log-file=05_cr2res_util_extract.log --output-dir=05_cr2res_util_extract_out cr2res_util_extract --smooth_slit=3 -smooth_spec=2.0E-7 05_cr2res_util_extract.sof 
echo "Next ${SETT}: cr2res_util_normflat 06_cr2res_util_normflat.sof"
esorex --log-file=06_cr2res_util_normflat.log --output-dir=06_cr2res_util_normflat_out cr2res_util_normflat 06_cr2res_util_normflat.sof

echo "Plotting trace and slit tilt"
show_trace.py 03_cr2res_util_trace_out/cr2res_util_calib_calibrated_collapsed_tw.fits 02_cr2res_util_calib_out/cr2res_util_calib_calibrated_collapsed.fits
plotsolution.py 03_cr2res_util_trace_out/cr2res_util_calib_calibrated_collapsed_tw.fits
show_raw.py -10,30 01_cr2res_cal_dark_out/cr2res_cal_dark_*_master.fits
show_raw.py 0.85,1.15 06_cr2res_util_normflat_out/cr2res_util_normflat_Open_master_flat.fits
if [ -f 04_cr2res_util_slit_curv_out/cr2res_util_calib_calibrated_collapsed_tw_tw.fits ] ; then
    show_trace_curv.py 04_cr2res_util_slit_curv_out/cr2res_util_calib_calibrated_collapsed_tw_tw.fits `grep WAVE 04_cr2res_util_slit_curv.sof | awk '{print $1}'`
fi


if [[ $SETT == K* ]] ||  [[ $SETT == H* ]] ; then
    CATNAME="redman"
else
    CATNAME="sarmiento"
fi
echo ${CATNAME}


if [[ $SETT == L* ]] ||  [[ $SETT == M* ]] ; then
    : # do nothing
else
    echo "Next ${SETT}: cr2res_util_calib 07_cr2res_util_calib.sof (calib+combine UNe raw frames)"
    esorex --log-file=07_cr2res_util_calib.log --output-dir=07_cr2res_util_calib_out cr2res_util_calib --collapse="MEAN" --subtract_nolight_rows=TRUE 07_cr2res_util_calib.sof
    echo "Next ${SETT}: cr2res_util_extract 08_cr2res_util_extract.sof (extract UNe)"
    esorex --log-file=08_cr2res_util_extract.log --output-dir=08_cr2res_util_extract_out cr2res_util_extract --smooth_slit=3 --swath_width=800 --oversample=4 08_cr2res_util_extract.sof
    echo "Next ${SETT}: cr2res_util_genlines 09_cr2res_util_genlines.sof"
    esorex --log-file=09_cr2res_util_genlines.log  --output-dir=09_cr2res_util_genlines_out  cr2res_util_genlines 09_cr2res_util_genlines.sof

    echo "Next ${SETT}: cr2res_util_wave 11_cr2res_util_wave.sof"
    esorex --log-file=11_cr2res_util_wave.log --output-dir=11_cr2res_util_wave_out cr2res_util_wave --wl_method=XCORR --wl_degree=0 --keep --wl_err=0.1 --fallback 11_cr2res_util_wave.sof
    #rm -f 11_cr2res_util_wave_out/cr2res_util_calib_calibrated_collapsed_extr1D_tw.fits
    #cp /home/tom/pCOMM/manucal/${SETT}_manucal_tw.fits 11_cr2res_util_wave_out/cr2res_util_calib_calibrated_collapsed_extr1D_tw.fits

    echo "Next ${SETT}: cr2res_util_wave 12_cr2res_util_wave.sof"
    esorex --log-file=12_cr2res_util_wave.log --output-dir=12_cr2res_util_wave_out cr2res_util_wave --wl_method=XCORR --wl_degree=2 --wl_err=0.03 --fallback 12_cr2res_util_wave.sof
    echo "Next ${SETT}: cr2res_util_calib 13_cr2res_util_calib.sof (calib+combine Etalon raw frames)"
    esorex --log-file=13_cr2res_util_calib.log --output-dir=13_cr2res_util_calib_out cr2res_util_calib --collapse="MEAN" 13_cr2res_util_calib.sof
    echo "Next ${SETT}: cr2res_util_extract 14_cr2res_util_extract.sof (extract Etalon)"
    esorex --log-file=14_cr2res_util_extract.log --output-dir=14_cr2res_util_extract_out cr2res_util_extract --swath_width=800 --oversample=4 --smooth_slit=3 14_cr2res_util_extract.sof
    echo "Next ${SETT}: cr2res_util_wave 15_cr2res_util_wave.sof"
    esorex --log-file=15_cr2res_util_wave.log --output-dir=15_cr2res_util_wave_out cr2res_util_wave --wl_method=ETALON --wl_degree=4 --fallback 15_cr2res_util_wave.sof

    show_wavecal.py 08_cr2res_util_extract_out/cr2res_util_calib_calibrated_collapsed_extr1D.fits 09_cr2res_util_genlines_out/lines_u_${CATNAME}.fits 09_cr2res_util_genlines_out/lines_u_${CATNAME}_une_sel.fits
    show_raw.py 07_cr2res_util_calib_out/cr2res_util_calib_calibrated_collapsed.fits # UNe
    show_raw.py 13_cr2res_util_calib_out/cr2res_util_calib_calibrated_collapsed.fits # Etalon

    for wavstep in 11 12 15 ; do
        plotsolution.py ${wavstep}_cr2res_util_wave_out/cr2res_util_calib_calibrated_collapsed_extr1D_tw.fits
    done
    plotresiduals.py -5e-3,5e-3 15_cr2res_util_wave_out/cr2res_util_calib_calibrated_collapsed_extr1D_lines_diagnostics.fits
    for wavstep in 11 12 ; do
        show_wavecal.py ${wavstep}_cr2res_util_wave_out/cr2res_util_calib_calibrated_collapsed_extr1D_extracted.fits 09_cr2res_util_genlines_out/lines_u_${CATNAME}.fits 09_cr2res_util_genlines_out/lines_u_${CATNAME}_une_sel.fits
    done
fi

exit 0;
