#! /opt/star/bin/perl -w

use File::Basename;
use Getopt::Std;
use Cwd 'abs_path';     # aka realpath()
my $resubmit = "resub.sh";

@jobDirs = (
"/star/u/brian40/PicoDst/PicoDst_NoMaker/test_KFParticle_gamma/condor_files"
);

@outputDirs = (
"/star/u/brian40/PicoDst/PicoDst_NoMaker/test_KFParticle_gamma/output/lam"
);

@flowCents = (
#"Data_160",
#"Data_161"
#,
#"Data_162"
#"Data_07",
#"Data_08",
#"Data_09_",
#"Data_10_",
#"Data_11_"
#"Data_12_",
#"Data_13_",
"Data"
#"Data_15_",
#"Data_16_"
);

@Centrality = (
"0",
"1",
"2",
"3",
"4",
"5",
"6",
"7",
"8"
);


foreach $eachCent (@flowCents) {
    
    foreach $outputDir (@outputDirs) {
    
        if (-e "$outputDir/$resubmit") { print `rm $outputDir/$resubmit \n`;}
    
        foreach $jobDir (@jobDirs) {
        
            $totalJobs=0.;
            $finishedJobs=0.;
            #$finishedJobs1=0.;
        
            foreach $centt (@Centrality) {
                
                foreach (glob ("$jobDir/${eachCent}/sched*.csh") ){
                    
                    my $eachScript = $_;
                    chomp $eachScript;
                    @fields = split(/\//,$eachScript) ;
                    $eachScriptNoPath= $fields[$#fields];
                    ($eachJOBID = $eachScriptNoPath) =~ s/\.csh// ;
                    chomp $eachJOBID;
                    
                    #if ( (!(-s "$outputDir/${eachCent}${centt}/$eachJOBID\_plam.root")) || (!(-s "$outputDir/${eachCent}${centt}/${eachJOBID}cen${centt}.weight\_112\_module\_new.root"))) {
                    #if ( (!(-s "$outputDir/${eachCent}${centt}/$eachJOBID\_plam.root")) || (!(-s "$outputDir/${eachCent}${centt}/${eachJOBID}cen${centt}.gamma112_fullEP\_eff\_pT02\_module.root"))) {
                    #if ( (!(-s "$outputDir/${eachJOBID}cen${centt}\_KFParticleQA.root")) || (!(-s "$outputDir/${eachJOBID}cen${centt}\_output.root")) || (!(-s "$outputDir/${eachJOBID}cen${centt}.v2_fullEP_eff_pT02_module.root")) || (!(-s "$outputDir/${eachJOBID}cen${centt}.weight_112_module_new.root"))) {
                    if (!(-s "$outputDir/${eachJOBID}cen${centt}.gamma112_fullEP\_eff\_pT02\_module.root")) { 
	 	    #if ( (!(-s "$outputDir/${eachJOBID}cen${centt}\_KFParticleQA.root")) || (!(-s "$outputDir/${eachJOBID}cen${centt}\_output.root")) || (!(-s "$outputDir/${eachJOBID}\_tree.root"))){	
			#print "$outputDir/${eachCent}${centt}/${eachJOBID}cen0.weight\_112\_module\_new.root";
                        
                        if (-e "$jobDir/${eachCent}/$eachJOBID.run.log") { print `rm $jobDir/${eachCent}/$eachJOBID.run.log \n`; }
                        
                        print "condor_submit $jobDir\/${eachCent}\/$eachJOBID.condor \n";
                        $submitCommand = "cd $jobDir\/${eachCent}; condor_submit $jobDir\/${eachCent}\/$eachJOBID.condor";
	 		$submitCommand1 = "sleep 3s";                      
 
                        open(macroFile,">>$outputDir/$resubmit");
                        print macroFile "$submitCommand \n";
			print macroFile "$submitCommand1 \n";
                    }
                }
                
                my @tempScripts = glob( "$jobDir/${eachCent}/sched*.csh" );
                #my @tempFinished = glob( "$outputDir/*KFParticleQA.root" ) ;
                #my @tempFinished1 = glob( "$outputDir/*output.root" ) ;
                #my @tempFinished2 = glob( "$outputDir/*v2_fullEP_eff_pT02_module.root" ) ;
                #my @tempFinished3 = glob( "$outputDir/*weight_112_module_new.root" ) ;
		my @tempFinished = glob( "$outputDir/*gamma112_fullEP\_eff\_pT02\_module.root");                
		#my @tempFinished2 = glob( "$outputDir/*_tree.root");

                #print scalar(@tempFinished1);
                
                $totalJobs +=scalar(@tempScripts);
                $finishedJobs +=scalar(@tempFinished);
                #$finishedJobs1 +=scalar(@tempFinished1);
                #$finishedJobs2 +=scalar(@tempFinished2);
                #$finishedJobs3 +=scalar(@tempFinished3);
            }
            
            my $perCentFinished = 0.;
            #my $perCentFinished1 = 0.;
            #my $perCentFinished2 = 0.;
            #my $perCentFinished3 = 0.;
            if ($totalJobs != 0) {$perCentFinished = ($finishedJobs/$totalJobs)*100.;}
            #if ($totalJobs != 0) {$perCentFinished1 = ($finishedJobs1/$totalJobs)*100.;}
            #if ($totalJobs != 0) {$perCentFinished2 = ($finishedJobs2/$totalJobs)*100.;}
            #if ($totalJobs != 0) {$perCentFinished3 = ($finishedJobs3/$totalJobs)*100.;}
            print "for $jobDir/${eachCent}, finished $finishedJobs jobs out of $totalJobs, ".$perCentFinished."% completed \n";
            #print "for $jobDir/${eachCent}, finished $finishedJobs1 jobs out of $totalJobs, ".$perCentFinished1."% completed \n";
            #print "for $jobDir/${eachCent}, finished $finishedJobs2 jobs out of $totalJobs, ".$perCentFinished2."% completed \n";
            #print "for $jobDir/${eachCent}, finished $finishedJobs3 jobs out of $totalJobs, ".$perCentFinished3."% completed \n";
        }
    }
}

exit;

