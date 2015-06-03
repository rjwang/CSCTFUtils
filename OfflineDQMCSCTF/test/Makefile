all:
	@echo "make run/fit"
run:
	@rm csctf.log
	@cmsRun run_offline_dqm_csctf_cfg.py > csctf.log
fit:
	@root -l -b -q readLCTAnalyzer.cc++ > output_lct.log
plot:
	@root -l -b -q plotter.cc
	@rm ~/www/my_public/CSCTF_LUTRel/*.png
	@rm ~/www/my_public/CSCTF_LUTRel/*.pdf
	@mv *.png ~/www/my_public/CSCTF_LUTRel/
	@mv *.pdf ~/www/my_public/CSCTF_LUTRel/
