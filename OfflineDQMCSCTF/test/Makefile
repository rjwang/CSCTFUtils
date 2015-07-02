all:
	@echo "make run/fit"
run:
	@cmsRun run_offline_dqm_csctf_cfg.py > /tmp/rewang/csctf.log &
test:
	@#cmsRun run_offline_dqm_csctf_cfg.py
	@root -l -b -q readLCTAnalyzer_MC_OneHitperStation.cc'("Input_MC_v1.root","output_test_1.root",0)'++
fit:
	@#root -l -b -q readLCTAnalyzer.cc'("Data13TeV_Run2015A_RAW.root","output_test.root")'++ > /tmp/rewang/output_lct.log
	@#root -l -b -q readLCTAnalyzer.cc'("Data13TeV_Run2015A_SingleMu.root","output_test.root")'++ > /tmp/rewang/output_lct.log
	@#root -l -b -q readLCTAnalyzer.cc'("swapCFEB1_2_rear_ME11b.root","output_swapCFEB1_2_rear_ME11b.root")'++ > /tmp/rewang/output_lct_swapCFEB1_2_rear_ME11b.log &
	@#root -l -b -q readLCTAnalyzer.cc'("swapCFEB3_4_rear_ME11b.root","output_swapCFEB3_4_rear_ME11b.root")'++ > /tmp/rewang/output_lct_swapCFEB3_4_rear_ME11b.log &
	@#root -l -b -q readLCTAnalyzer.cc'("swapCFEB1_2_front_ME11b.root","output_swapCFEB1_2_front_ME11b.root")'++ > /tmp/rewang/output_lct_swapCFEB1_2_front_ME11b.log &
	@#root -l -b -q readLCTAnalyzer.cc'("swapCFEB3_4_front_ME11b.root","output_swapCFEB3_4_front_ME11b.root")'++ > /tmp/rewang/output_lct_swapCFEB3_4_front_ME11b.log &
	@#root -l -b -q readLCTAnalyzer.cc'("/tmp/rewang/csctf_Run246963_ME-11_match_ME+11.root","output_test.root")'++ > /tmp/rewang/output_lct.log
	@#root -l -b -q readLCTAnalyzer.cc'("/tmp/rewang/csctf_Run246963_nochange.root","output_test.root")'++ > /tmp/rewang/output_lct.log
	@root -l -b -q readLCTAnalyzer.cc'("raw_test.root","output_test.root")'++ > /tmp/rewang/output_lct.log &

plot:
	@hadd -f output_test.root output_test_*.root
	@root -l -b -q plotter.cc'("output_test.root")'
	@rm ~/www/my_public/CSCTF_LUTRel/*.png
	@rm ~/www/my_public/CSCTF_LUTRel/*.pdf
	@mv *.png ~/www/my_public/CSCTF_LUTRel/
	@mv *.pdf ~/www/my_public/CSCTF_LUTRel/

clean:
	@rm output_test_*.root
