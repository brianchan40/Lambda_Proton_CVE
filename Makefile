#CENTRALITIES = 5 6 7 8
#CENTRALITIES = 4
CENTRALITIES = 0 1 2
#PERIOD_5 = 160 161 162
#PERIOD_5 = _low_ _mid_ _high_ _
#PERIOD_5 = 10
PERIOD_5 = 14
#PERIOD_5 = 08 09 10 11 12 13 14 15 16 17
PERIOD_4 = 1 2 3 4 5 6 7 8 9
#PERIOD_4 = 1
#LAMBDA = lam antilam
LAMBDA = lam
#LUM = high nor mid low
LUM = low

submit_analysis:
	rm -rf ./pla*.zip  ./pla*.package
	rm -rf ./output/*
	for lam in $(LAMBDA) ; do ./submit_scheduler.sh $$lam 0 8 ; done ;

submit_check:
	rm -rf ./plac_*.zip  ./plac_*.package
	rm -rf ./output/combine/
	for cen in $(CENTRALITIES) ; do for beginnum in $(BEGINNUM) ; do ./submit_combine.sh 09 $$cen $$beginnum ; done ; done ;

c:
	rm ./condor_files/scripts/*0.csh 
	rm ./condor_files/scripts/*1.csh 
	rm ./condor_files/scripts/*2.csh
	rm ./condor_files/scripts/*3.csh
	rm ./condor_files/scripts/*4.csh
	rm ./condor_files/scripts/*5.csh
	rm ./condor_files/scripts/*6.csh
	rm ./condor_files/scripts/*7.csh 
	rm ./condor_files/scripts/*8.csh
	rm ./condor_files/scripts/*9.csh
	rm ./condor_files/*0.list
	rm ./condor_files/*1.list
	rm ./condor_files/*2.list
	rm ./condor_files/*3.list
	rm ./condor_files/*4.list
	rm ./condor_files/*5.list
	rm ./condor_files/*6.list
	rm ./condor_files/*7.list
	rm ./condor_files/*8.list
	rm ./condor_files/*9.list
	rm -rf condor_files/*
	mkdir condor_files/scripts
	rm sched*
	rm -r output/*_debug

cc:
	rm -rf condor_files_check/*
	rm -r output_check/*_debug

combine:
	for lam in $(LAMBDA) ; do for cen in $(CENTRALITIES) ; do \
                hadd ./Results_$${lam}_18/cen$$cen.weight_112_module.root ./temp/*/*cen$$cen.weight_112_module_new.root; \
                hadd ./Results_$${lam}_18/cen$$cen.v2_fullEP_eff_pT02_module.root ./temp/*/*cen$$cen.v2_fullEP_eff_pT02_module.root; \
		hadd ./Results_$${lam}_18/cen$${cen}_plam.root ./temp/*/*cen$${cen}_output.root; \
        done; done;

combine1:
	for lam in $(LAMBDA) ; do for cen in $(CENTRALITIES) ; do \
                mv ./Results_$${lam}_18/cen$$cen.weight_112_module.root ./Results_$${lam}_18/cen$$cen.weight_112_module_old.root; \
                mv ./Results_$${lam}_18/cen$$cen.v2_fullEP_eff_pT02_module.root ./Results_$${lam}_18/cen$$cen.v2_fullEP_eff_pT02_module_old.root; \
                mv ./Results_$${lam}_18/cen$${cen}_plam.root ./Results_$${lam}_18/cen$${cen}_plam_old.root; \
                hadd ./Results_$${lam}_18/cen$$cen.weight_112_module.root ./temp/*/*cen$$cen.weight_112_module_new.root; \
                hadd ./Results_$${lam}_18/cen$$cen.v2_fullEP_eff_pT02_module.root ./temp/*/*cen$$cen.v2_fullEP_eff_pT02_module.root; \
                hadd ./Results_$${lam}_18/cen$${cen}_plam.root ./temp/*/*cen$${cen}_output.root; \
        done; done ;

combine2:
	for lam in $(LAMBDA) ; do for lum in $(LUM) ; do for cen in $(CENTRALITIES) ; do \
        	mv ./Results_$${lam}_18/cen$${cen}_plam.root ./Results_$${lam}_18/cen$${cen}_plam_old2.root; \
                hadd ./Results_$${lam}_18/cen$${cen}_plam.root ./temp/*/*cen$${cen}_plam.root; \
                hadd ./Results_$${lam}_18/cen$$cen.gamma112_fullEP_eff_pT02_module.root ./temp/*/*cen$$cen.gamma112_fullEP_eff_pT02_module.root; \
	done; done; done;

combine3:
	for lam in $(LAMBDA) ; do for lum in $(LUM) ; do for cen in $(CENTRALITIES) ; do \
		hadd ./Results_$${lam}_$${lum}_14/cen$$cen.gamma112_fullEP_eff_pT02_module.root ./temp/*cen$$cen.gamma112_fullEP_eff_pT02_module.root; \
                hadd ./Results_$${lam}_$${lum}_14/cen$$cen.plam_final.root ./temp/*cen$${cen}plam.root; \
        done; done; done;

combine_2:
	for lam in $(LAMBDA) ; do for lum in $(LUM) ; do \
		hadd 2014_$${lum}_info_$${lam}.root ./output_check/$${lam}_$${lum}/*/*2014_high_info_1.root; \
	done ; done ;


prep:
	for daynum in $(PERIOD_5) ; do for cen in $(CENTRALITIES) ; do mkdir ./condor_files/Data ; done; done ;

dataquery:
	get_file_list.pl -keys 'path,filename' -cond 'production=P17ih,trgsetupname=AuAu54_production_2017,filetype=daq_reco_MuDst,filename~st_physics_181,storage=HPSS' -limit 500
kill:
	condor_rm brian40	
check:
	@condor_q brian40 | tail -10
test:
	rm -rf ./pla_*.zip  ./pla_*.package
	./submit_scheduler.sh 161 3 gamma 300 0.1
	./submit_scheduler.sh 156 3 gamma 300 0.1
test_june:
	rm -rf ./pla_*.zip  ./pla_*.package
	./submit_scheduler.sh 162 3 gamma 800 0.1
	#./submit_scheduler.sh 161 3 gamma 500 0.1
	#./submit_scheduler.sh 163 3 gamma 300 0.1
test_phi:
	rm -rf ./pla_*.zip  ./pla_*.package
	./submit_scheduler.sh 161 3 phi 300 0.1
	./submit_scheduler.sh 156 3 phi 300 0.1
test_loop:

.PHONY: clean submit_phi submit_ep dataquery kill check test test_loop
