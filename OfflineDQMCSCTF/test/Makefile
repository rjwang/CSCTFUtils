all:
	@echo "make run/fit"
run:
	@rm csctf.log
	@cmsRun run_offline_dqm_csctf_cfg.py > /tmp/rewang/csctf.log
test:
	@cmsRun run_offline_dqm_csctf_cfg.py
fit:
	@#root -l -b -q readLCTAnalyzer.cc'("Data13TeV_Run2015A_RAW.root","output_test.root")'++ > /tmp/rewang/output_lct.log
	@#root -l -b -q readLCTAnalyzer.cc'("Data13TeV_Run2015A_SingleMu.root","output_test.root")'++ > /tmp/rewang/output_lct.log
#	@root -l -b -q readLCTAnalyzer.cc'("/tmp/rewang/csctf_Run246963_ME-11_match_ME+11.root","output_test.root")'++ > /tmp/rewang/output_lct.log
	@root -l -b -q readLCTAnalyzer.cc'("/tmp/rewang/csctf_Run246963_nochange.root","output_test.root")'++ > /tmp/rewang/output_lct.log



plot:
	@root -l -b -q plotter.cc
	@rm ~/www/my_public/CSCTF_LUTRel/*.png
	@rm ~/www/my_public/CSCTF_LUTRel/*.pdf
	@mv *.png ~/www/my_public/CSCTF_LUTRel/
	@mv *.pdf ~/www/my_public/CSCTF_LUTRel/
