#include <ctime>

#include "BatchScriptWriter.h"

BatchScriptWriter::BatchScriptWriter(int argc, char* argv[]):
  exec("")
{
  for (int i=0; i<argc; i++){
    // skip arguments we dont need
    if ( string(argv[i])==string("--nbatchjobs") || (i>0 && string(argv[i-1])==string("--nbatchjobs")) ||
         string(argv[i])==string("--batchstartn") || (i>0 && string(argv[i-1])==string("--batchstartn")) ||
         string(argv[i])==string("--batcheos")
           ) {
      continue;
    }
    exec += string(argv[i]) + " ";
  }
  subpkg = string(argv[0]);
}

BatchScriptWriter::~BatchScriptWriter(){}

void BatchScriptWriter::writeScripts(OptParser *arg, vector<Combiner*> *cmb){

  for ( int i=0; i<arg->combid.size(); i++ )
  {
    int combinerId = arg->combid[i];
    Combiner *c = cmb->at(combinerId);
    c->combine();
    if ( !c->isCombined() ) continue;

    TString methodname = "";
    if ( arg->isAction("pluginbatch") ) {
      methodname = "Plugin";
    }
    else if ( arg->isAction("coveragebatch") ) {
      methodname = "Coverage";
    }
    else {
      cout << "BatchScriptWriter::writeScripts() : ERROR : only need to write batch scripts for pluginbatch method" << endl;
      exit(1);
    }
    if ( arg->isAction("uniform") ) methodname += "Uniform";
    if ( arg->isAction("gaus") ) methodname += "Gaus";

    cout << "Writing submission scripts for combination " << c->getName() << endl;

    TString dirname = "scan1d"+methodname+"_"+c->getName()+"_"+arg->var[0];
    if ( arg->var.size()==2 ) {
      dirname = "scan2d"+methodname+"_"+c->getName()+"_"+arg->var[0];
    }
    if ( arg->var.size()>1) {
      dirname += "_"+arg->var[1];
    }
    if ( arg->isAction("coveragebatch") ) {
      dirname += arg->id<0 ? "_id0" : Form("_id%d",arg->id);
    }
    TString scripts_dir_path = "sub/" + dirname;
    TString outf_dir = "root/" + dirname;
    system(Form("mkdir -p %s",scripts_dir_path.Data()));
    TString scriptname = "scan1d"+methodname+"_"+c->getName()+"_"+arg->var[0];
    if ( arg->isAction("coveragebatch") ) {
      scriptname += arg->id<0 ? "_id0" : Form("_id%d",arg->id) ;
    }
    // if write to eos then make the directory
    if ( arg->batcheos ) {
      time_t t = time(0);
      struct tm * now = localtime(&t);
      int day = now->tm_mday;
      int month = now->tm_mon+1;
      int year  = now->tm_year+1900;

      std::string name(getenv("USER"));
      char initial = name[0];
      std::cout << "Name: " << name << " first: " << initial <<std::endl;
      TString eos_path = Form("/eos/lhcb/user/%c/%s/gammacombo/%02d%02d%04d",initial,name.c_str(),day,month,year);
      // TString eos_path = Form("/eos/lhcb/user/t/tmombach/gammacombo/%02d%02d%04d",day,month,year);
      // system(Form("/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select mkdir %s",eos_path.Data()));
      eos_path += Form("/%s",dirname.Data());
      // system(Form("/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select mkdir -p %s",eos_path.Data()));
      outf_dir = eos_path;
    }
    if ( arg->var.size()==2 ) {
      scriptname = "scan2d"+methodname+"_"+c->getName()+"_"+arg->var[0];
    }
    if ( arg->var.size()>1) {
      scriptname += "_"+arg->var[1];
    }
    scriptname = scripts_dir_path + "/" + scriptname;

    // write out a list of the submission files
    ofstream subfilelist;
    subfilelist.open(scriptname + "_sublist.txt");
    for ( int job=arg->batchstartn; job<arg->batchstartn+arg->nbatchjobs; job++ ) {
      TString fname = scriptname + Form("_run%d",job) + ".sh";
      subfilelist << fname << endl;
      writeScript(fname, outf_dir, job, arg);
    }
    subfilelist.close();
    cout << "Written submission file list to\n\t" << scriptname << "_sublist.txt" << endl;
    writeCondorScript(scriptname, arg);
  }
}

// Copy the same submission script with minor differences for the datasets option
// -> only difference is that the NLL is stored in the pdf already and does not have to be combined
// -> no combiners for the datasets option
void BatchScriptWriter::writeScripts_datasets(OptParser *arg, PDF_Abs* pdf){

  if(!pdf){
    std::cout<< "BatchScriptWriter::writeScripts_datasets(): ERROR: No PDF given " << std::endl;
    exit(1);
  }

  TString methodname = "";
  if ( arg->isAction("pluginbatch") ) {
    methodname = "Plugin";
  }
  else if ( arg->isAction("coveragebatch") ) {
    methodname = "Coverage";
  }
  else {
    cout << "BatchScriptWriter::writeScripts() : ERROR : only need to write batch scripts for pluginbatch method" << endl;
    exit(1);
  }
  if ( arg->isAction("uniform") ) methodname += "Uniform";
  if ( arg->isAction("gaus") ) methodname += "Gaus";

  cout << "Writing submission scripts for PDF " << pdf->getName() << endl;

  TString dirname = "scan1dDatasets"+methodname+"_"+pdf->getName()+"_"+arg->var[0];
  if ( arg->var.size()==2 ) {
    dirname = "scan2dDatasets"+methodname+"_"+pdf->getName()+"_"+arg->var[0];
  }
  if ( arg->var.size()>1) {
    dirname += "_"+arg->var[1];
  }
  if ( arg->isAction("coveragebatch") ) {
    dirname += arg->id<0 ? "_id0" : Form("_id%d",arg->id);
  }
  TString scripts_dir_path = "sub/" + dirname;
  TString outf_dir = "root/" + dirname;
  system(Form("mkdir -p %s",scripts_dir_path.Data()));
  TString scriptname = "scan1dDatasets"+methodname+"_"+pdf->getName()+"_"+arg->var[0];
  if ( arg->isAction("coveragebatch") ) {
    scriptname += arg->id<0 ? "_id0" : Form("_id%d",arg->id) ;
  }
  // if write to eos then make the directory
  if ( arg->batcheos ) {
    time_t t = time(0);
    struct tm * now = localtime(&t);
    int day = now->tm_mday;
    int month = now->tm_mon+1;
    int year  = now->tm_year+1900;

    std::string name(getenv("USER"));
    char initial = name[0];
    std::cout << "Name: " << name << " first: " << initial <<std::endl;
    TString eos_path = Form("/eos/lhcb/user/%c/%s/gammacombo/%02d%02d%04d",initial,name.c_str(),day,month,year);
    // TString eos_path = Form("/eos/lhcb/user/t/tmombach/gammacombo/%02d%02d%04d",day,month,year);
    //system(Form("/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select mkdir %s",eos_path.Data()));
    // system(Form("mkdir %s",eos_path.Data()));
    eos_path += Form("/%s",dirname.Data());
    //system(Form("/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select mkdir -p %s",eos_path.Data()));
    // system(Form("mkdir -p %s",eos_path.Data()));
    outf_dir = eos_path;
  }
  if ( arg->var.size()==2 ) {
    scriptname = "scan2dDatasets"+methodname+"_"+pdf->getName()+"_"+arg->var[0];
  }
  if ( arg->var.size()>1) {
    scriptname += "_"+arg->var[1];
  }
  scriptname = scripts_dir_path + "/" + scriptname;

  // write out a list of the submission files
  ofstream subfilelist;
  subfilelist.open(scriptname + "_sublist.txt");
  for ( int job=arg->batchstartn; job<arg->batchstartn+arg->nbatchjobs; job++ ) {
    TString fname = scriptname + Form("_run%d",job) + ".sh";
    subfilelist << fname << endl;
    writeScript(fname, outf_dir, job, arg);
  }
  subfilelist.close();
  writeCondorScript(scriptname, arg);
}

void BatchScriptWriter::writeCondorScript(TString fname, OptParser *arg) {

  TString subfilename = fname + ".sub";
  //cout << "\t" << subfilename << endl;
  ofstream outfile;
  outfile.open(subfilename);

  char cwd[1024];
  getcwd(cwd,1024);

  outfile << "##### auto-generated by BatchScriptWriter #####" << endl;
  outfile << "executable = $(subfile)" << endl;
  outfile << "arguments = " << endl;
  //outfile << "requirements = OpSysAndVer == \"CentOS7\"" << endl;
  outfile << "output = $(subfile).out" << endl;
  outfile << "error = $(subfile).err" << endl;
  outfile << Form("log = %s.log",fname.Data()) << endl;
  if (arg->queue != "") outfile << "+JobFlavour = " << '"' << arg->queue << '"'<< endl;
  outfile << Form("queue subfile from %s_sublist.txt",fname.Data()) << endl;
  outfile.close();

  system(Form("chmod +x %s",subfilename.Data()));

  cout << "Written condor submission script to\n\t" << subfilename << endl;

  if ( arg->queue != "" ) {
    //cout << Form("condor_submit %s/%s",cwd,subfilename.Data()) << endl;
    system(Form("condor_submit %s/%s",cwd,subfilename.Data()));
  }
}


void BatchScriptWriter::writeScript(TString fname, TString outfloc, int jobn, OptParser *arg) {

  TString rootfilename = fname;
  (rootfilename.ReplaceAll("sub","root")).ReplaceAll(".sh",".root");
  cout << "\t" << fname << endl;
  ofstream outfile;
  outfile.open(fname);

  char cwd[1024];
  getcwd(cwd,1024);

  outfile << "#!/bin/bash" << endl;
  outfile << "##### auto-generated by BatchScriptWriter #####" << endl;
  outfile << Form("rm -f %s/%s.done",cwd,fname.Data()) << endl;
  outfile << Form("rm -f %s/%s.fail",cwd,fname.Data()) << endl;
  outfile << Form("rm -f %s/%s.run",cwd,fname.Data()) << endl;
  outfile << Form("rm -f %s/%s.log",cwd,fname.Data()) << endl;
  outfile << "mkdir -p scratch" << endl;
  outfile << "cd scratch" << endl;
  outfile << Form("source %s/../scripts/setup_lxplus.sh",cwd) << endl;
  outfile << Form("cp -r %s/ExpNll .",cwd) << endl;
  outfile << "mkdir -p bin" << endl;
  outfile << Form("cp %s/%s bin/",cwd,subpkg.c_str()) << endl;
  outfile << "mkdir -p plots/dot" << endl;
  outfile << Form("cp -r %s/plots/dot/* plots/dot",cwd) << endl;
  outfile << "mkdir -p plots/par" << endl;
  outfile << Form("cp -r %s/plots/par/* plots/par",cwd) << endl;
  outfile << "mkdir -p plots/scanner" << endl;
  outfile << "mkdir -p root" << endl;
  outfile << Form("touch %s/%s.run",cwd,fname.Data()) << endl;
  outfile << Form("if ( %s --nrun %d ); then",exec.c_str(),jobn) << endl;
  outfile << Form("\ttouch %s/%s.done",cwd,fname.Data()) << endl;
  outfile << Form("\trm -f %s/%s.run",cwd,fname.Data()) << endl;
  outfile << "else" << endl;
  outfile << Form("\ttouch %s/%s.fail",cwd,fname.Data()) << endl;
  outfile << Form("\trm -f %s/%s.run",cwd,fname.Data()) << endl;
  outfile << "fi" << endl;
  
  TString basename = rootfilename;
  basename.Remove(0, rootfilename.Last('/')+1);
  system(Form("mkdir -p %s/%s",cwd,outfloc.Data()));
  TString copy_line = Form("cp %s %s/%s/%s",rootfilename.Data(),cwd,outfloc.Data(),basename.Data());
  if ( arg->batcheos ) copy_line = Form("xrdcp -p %s root://eoslhcb.cern.ch/%s/%s",rootfilename.Data(),outfloc.Data(),basename.Data());
  outfile << copy_line << endl;

  outfile.close();

  system(Form("chmod +x %s",fname.Data()));

  // if ( arg->queue != "" ) {
  //   //system(Form("bsub -R \"rusage[mem=40000]\" -q %s -o %s/%s.log %s/%s",arg->queue.Data(),cwd,fname.Data(),cwd,fname.Data()));
  //   system(Form("bsub -q %s -o %s/%s.log %s/%s",arg->queue.Data(),cwd,fname.Data(),cwd,fname.Data()));
  // }
}
