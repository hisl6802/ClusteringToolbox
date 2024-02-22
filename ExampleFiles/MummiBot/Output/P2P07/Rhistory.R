# PID of current job: 2469780
mSet<-InitDataObjects("mass_all", "mummichog", FALSE)
mSet<-SetPeakFormat(mSet, "rmp")
mSet<-UpdateInstrumentParameters(mSet, 5.0, "positive", "yes", 0.02);
mSet<-Read.PeakListData(mSet, "Replacing_with_your_file_path");
mSet<-SanityCheckMummichogData(mSet)
mSet<-SetPeakEnrichMethod(mSet, "mum", "v2")
mSet<-SetMummichogPval(mSet, 0.05)
mSet<-PerformPSEA(mSet, "mmu_kegg", "current", 3 , 100)
mSet<-PlotPeaks2Paths(mSet, "peaks_to_paths_8_", "png", 72, width=NA)
mSet<-SaveTransformedData(mSet)
