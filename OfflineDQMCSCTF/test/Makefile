all:
	@echo "make run/fit"
run:
	@rm csctf.log
	@cmsRun run_offline_dqm_csctf_cfg.py > /tmp/rewang/csctf.log
test:
	@cmsRun run_offline_dqm_csctf_cfg.py
fit:
	@#root -l -b -q readLCTAnalyzer.cc'("csctf_Run246963.root","output_test.root",true)'++ > /tmp/rewang/output_lct.log ## plus endcap
	@root -l -b -q readLCTAnalyzer.cc'("csctf_Run246963.root","output_test.root",false)'++ > /tmp/rewang/output_lct.log ## minus endcap
plot:
	@root -l -b -q plotter.cc
	@rm ~/www/my_public/CSCTF_LUTRel/*.png
	@rm ~/www/my_public/CSCTF_LUTRel/*.pdf
	@mv *.png ~/www/my_public/CSCTF_LUTRel/
	@mv *.pdf ~/www/my_public/CSCTF_LUTRel/
