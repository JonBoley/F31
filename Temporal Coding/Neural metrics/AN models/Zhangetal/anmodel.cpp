/** \file anmodel.cpp
 * \brief This is the main file for the AN model.
 * <p>the output is stored in filterout, synapseout
 * <p>The stimulus is read from the file(default filename is "click")
 * <p>The program initializes the filter bank model first(according to the parameters
 * specified in the commandline).
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include "cmpa.hpp"
#include "complex.hpp"

#define NSTIMMAX 100000 // enough points for durations up to 250 msec
#define NUMCHANS 64 // Number of fibers in the population to be simulated
#define MAXSPIKES 10000
double stim[NSTIMMAX];
double sout[NSTIMMAX];
double sptime[MAXSPIKES]; //spike time storage

/// This structure defines the parameters needed for the model
typedef struct {
  char wavefile[100]; //filename for stimulus waveform
  int banks; //number of filters in filterbank
  double delx;  //"delta x", or spacing of "IHC"s along BM (linearly spaced IHCs in terms of location along BM
				//  results in approximately logarithmic spacing in terms of frequency.
  int stim;
  double tdres; //time step (time domain resolution)
  int nstim,nrep; // nstim=number of pts in stimulus; nrep is number of reps of stimulus
  int runmode;
  double freqt,spl,dur,reptim;    //for tone and other stim
                                  //noise(2),click(3)
  double F1,F2,spl1,spl2;         //for two tone suppression(8)
  double *buf;
  double rf;                      //rise & fall
  int rfmethod;// rise & fall method
  double cfhi,cflo,cfinc;
  double splhi,spllo,splinc;
  double freqthi,freqtlo,freqtinc;
} T_stim;

int parsecommandline(TAuditoryNerve *pan,T_stim *ptm,int argc, char* argv[]);
int synapse(TAuditoryNerve *pan,int nstim,int nrep,double *sout,double *sptime);

int main(int argc, char* argv[])
{
  static T_stim tm;                    // structure to store the input parameters
  static TAuditoryNerve an[NUMCHANS];              // filter bank(Maximum 64)
  double xlo,xhi,delx,xcenter;

  int species = 1; //CAT
  int icenter;
  int icf;
  int nspikes;

  double x,cf;
  FILE *fpstim,*fpbm,*fpsout;          //open these file when run
  double buf; //tempory value to store the stimulus
  long flag,nowstim,totalstim;

  an[0].tdres = tm.tdres = 10e-6;      //default sample rate = 500kHz
  an[0].cf = 1000;                     //default CF      : 1000 Hz
  an[0].spont = 50;                    //     spont rate : 50
  strcpy(tm.wavefile,"click");         //input file      : click
  tm.reptim = 0.020;                   //repetition time  : 20 msec
  tm.delx = 0.05;                      //distance between two filters along BM : 0.05mm
  tm.cflo = -1;
  tm.cfhi = -1;
  tm.banks = -1;                       //default         : no filters
  tm.nrep = 0;

  parsecommandline(&an[0],&tm,argc,argv); // get the parameters from command line

  if(tm.banks == -1) //this is for a simulation of a single anditory nerve fiber
    {
      tm.banks = 1;
      xcenter = cochlea_f2x(species,an[0].cf); //convert from CF to BM location
      delx = 0;
    }
  else if(tm.cflo == -1) //this is for AN fiber population
    {
      xcenter = cochlea_f2x(species,an[0].cf);
      delx = tm.delx;
    }
  else
    {
      if(tm.cfhi<tm.cflo) error("cfhi should be greater than cflo");
      xlo = cochlea_f2x(species,tm.cflo);
      xhi = cochlea_f2x(species,tm.cfhi);
      xcenter = (xhi -xlo)/2. + xlo;
      delx = (xhi-xlo)/(tm.banks-1);
      an[0].cf = cochlea_x2f(species,xcenter);
    };

  printf("\n Filter Bank: center CF=%f, xcenter=%f(mm from apex), number of fibers= %d\n",an[0].cf,xcenter,tm.banks);

  icenter = tm.banks/2;

  nowstim = 0;
  totalstim = (long)(tm.reptim/an[0].tdres);
  buf = 0.0;

  // ******************now begin to pass the stimulus to the model************************
  // Read in the stimulus
  for(nowstim=0; nowstim<NSTIMMAX; nowstim++) stim[nowstim] = 0.0;
  printf("\nReading in stimulus \n");
  fpstim = fopen(tm.wavefile,"r");
  nowstim=0;
  while((fscanf(fpstim,"%le",&stim[nowstim]) != EOF) & nowstim<NSTIMMAX)
    nowstim++;
  printf("%d points read in\n",nowstim);
  if(totalstim < nowstim) totalstim=nowstim;
  fclose(fpstim);
  if(nowstim>=NSTIMMAX) error("Entire stim not read in: increase NSTIMMAX");
  if(totalstim>=NSTIMMAX) error("Entire reptim too long: increase NSTTMMAX");

  // construct the model
  for(icf = 0; icf < tm.banks; icf++) /* Loop through the filters */
    {
      x = xcenter + (icf - icenter)*delx;
      cf = cochlea_x2f(species,x); // from BM location to CF
      an[icf].cf = cf;
      an[icf].tdres = an[0].tdres;
      an[icf].spont = an[0].spont;
      an[icf].construct();
      an[icf].init(cf,an[icf].spont);
      printf("icf = %d out of %d x= %f cf= %f\n",icf,tm.banks,x,cf);
    };

  flag = 0; //if = 1,don't read from stimulus file but set zero as input
  nowstim = 0;
  fpbm = fopen("filterout","w");
  fpsout = fopen("synapseout","w");
  while(nowstim<totalstim)
    {
	  buf = stim[nowstim];
      fprintf(fpbm,"%g",nowstim*an[0].tdres);
      fprintf(fpsout,"%g",nowstim*an[0].tdres);
      for(icf = 0; icf < tm.banks; icf++) // Loop through the filters
		{
		  an[icf].run(buf);
		  fprintf(fpbm," %g",an[icf].bmbuf);  // Write filter response out to file
		  fprintf(fpsout," %g",an[icf].sbuf);  // Write synapse output to file
	  };
      fprintf(fpbm,"\n");
      fprintf(fpsout,"\n");
	  sout[nowstim] = an[0].sbuf;
      nowstim++;
    }; //end of while
  fclose(fpbm);
  fclose(fpsout);
  printf("\n AN filter bank simulation is complete. \n");
  if((tm.nrep>0)&&(tm.banks==1))
	{
	  printf("\n Generate spike times and the PST histogram(write out to files pstplot and sptime...\n");
      nspikes = synapse(&an[0],totalstim,tm.nrep,sout,sptime);
      double bin = 0.0001;
	  int nbins = (int)(tm.reptim / bin + 1);
	  int *pst = new int[nbins];
	  int i;
	  /* Construct PST HISTOGRAM, Interval hist, Period hist.  */
	  for(i = 0; i < nbins; i++) pst[i] = 0;
	  for(i = 0; i < nspikes; i++)
      {
      	/* both sptime and bin are in secs */
	    int ipst = (int)(fmod(sptime[i],tm.tdres*totalstim) / bin);
        pst[ipst] = pst[ipst] + 1;
      }
 	  fstream *fp = new fstream("pstplot",ios::out);
      for(i = 0; i < nbins; i++)
	  {
 		  	(*fp)<<(i*bin)<<' '<<pst[i]<<"\n";
      };
      delete fp;
      delete pst;

 	  fp = new fstream("sptime",ios::out);
      for(i = 0; i < nspikes; i++)
	  {
 		  	(*fp)<<i<<' '<<sptime[i]<<"\n";
      };
      delete fp;
	};
};

// ----------------------------------------------------------------
/**
  * This function reads the command line options and store them in the structure
*/
int parsecommandline(TAuditoryNerve *pan,T_stim *ptm,int argc, char* argv[])
{
  int i;
  int needhelp = 0;
  char *para;

  i = 1;
  if(argc ==1 ) needhelp = 1;
  while(i<argc)
  { para = argv[i];
    if((*para)=='-')
      { para++;
	if(strcmp(para,"tdres")==0)
	  {
	    pan->tdres = double(atof(argv[++i]));
	    ptm->tdres = pan->tdres;
	  }
	else if(strcmp(para,"cf")==0)
	  pan->cf = double(atof(argv[++i]));
	else if(strcmp(para,"spont")==0)
	  pan->spont = double(atof(argv[++i]));

	else if(strcmp(para,"fibers")==0)
	  ptm->banks = atoi(argv[++i]);
	else if(strcmp(para,"reptim")==0)
	  ptm->reptim = double(atof(argv[++i]));
	else if(strcmp(para,"delx")==0)
	  ptm->delx = double(atof(argv[++i]));
	else if(strcmp(para,"trials")==0)
	  ptm->nrep = atoi(argv[++i]);
	else if(strcmp(para,"cfhi")==0)
	  ptm->cfhi = double(atof(argv[++i]));
	else if(strcmp(para,"cflo")==0)
	  ptm->cflo = double(atof(argv[++i]));
	else if(strcmp(para,"wavefile")==0)
	  {
	    ptm->stim = 11;
	    strcpy(ptm->wavefile,argv[++i]);
	  }
	else if(strcmp(para,"help")==0)
	  needhelp = 1;
	else
	  { printf("\nUnkown parameters --> %s",para); needhelp = 1; break; };
      }
    else { printf("\nUnkown parameters --> %s",para); needhelp = 1; break; };
    i++;
  };
  if(needhelp==1)
    {
      printf("\n This program accepts following parameters:\n");

      printf("\n -cf #(1000)   --> character frequency of the AN model fiber (or center CF for filter bank)");
      printf("\n -spont #(50)  --> spontaneous rate of the fiber");
      printf("\n -tdres #(10e-6)--> time domain resolution(in seconds)");

      printf("\n -fibers #(1) --># of filters, use this option with cf,cflo,cfhi,delx");
      printf("\n -wavefile filename(click) --> specify the stimulus wavefile name(click)");
      printf("\n -delx #(0.05mm) --> distance between filters along basilar membrane");
      printf("\n -reptim #(0.02) --> Total duration of stimulation(20msec)");
      printf("\n -cfhi #(-1)   --> highest cf in population (specify cfhi,cflo will recalculate cf&delx");
      printf("\n -cflo #(-1)   --> lowest cf in population");
      printf("\n -trials #(0)  --> spike generation trials (this is valid only if # of fibers=1)");
      printf("\n");
      exit(1);
    };
return(0);
};

// --------------------------------------------------------------------------------
/** this routine generates spike times according to *sout with nrep representation
   put the spikes time into buffer *spikes, nspikes
*/
int synapse(TAuditoryNerve *pan,int nstim,int nrep,double *sout,double *sptime)
{
  int i,j,isp;
  int nspikes = 0;
  /* sout[] contains prob of spike vs. time - run through it nrep times to look for spikes */
  isp = 0;
  sptime[0] = 0.; /* start out with a spike at t=0 on 1st rep */
  pan->sg->init(sout[0]); //initialize the spike generator
  for(i = 0; i < nrep; i++)
    {
      pan->sg->init(sout[0]);
      for( j = 0; j < nstim; j++)
	{
	  if(pan->sg->run(sout[j])==1)
	    {
	      sptime[isp] = pan->sg->Get_rtime();  /* leave sptime in sec */
	      isp++;
		  if(isp>MAXSPIKES) error("\nToo much spikes generated: increase MAXSPIKES\n");
	    }
	}
      nspikes = isp;
    };
  return nspikes;
};
