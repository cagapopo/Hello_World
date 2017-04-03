#include "MyAna.hpp"

#include "PtOrder.hpp"

#include <TF1.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TKDE.h>
#include <TLatex.h>
#include <TEllipse.h>
#include <TLine.h>
#include <TProfile.h>
#include <assert.h>  
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <TRandom.h>
#include "DetectorGeometry.hpp"
#include "DetectorConstants.hpp"
#define PI 3.14159265359
using namespace std;






float CELERITY=(299.792458);//mm/ns
float TIMEOFFSET=11.8;

bool MAKETIMEPLOT=false;
bool EVTDISPLAY=false;
bool VBF = false;
bool MB = true;


int WIDTHWINDOW=20;//ps (see binning of h_time)
float ENERGYCUT=0.02;
float TIMECUT=1250;//12.5
float JETPT=30000;
float ETAMIN=2.5;
float ETAMAX=3.2;



float slidingWindow(TH1F*hist,int width)//width should be even
  {
    int theBin=0;
    float max=0;
    for(int ibin=width/2; ibin<hist->GetNbinsX()-1-width/2; ibin++)
      {
	float sum=0;
	for(int jbin=ibin-width/2;jbin<ibin+width/2;jbin++)
	  {
	    sum+=hist->GetBinContent(jbin);
	  }
	if(sum>max)
	  {
	    max=sum;
	    theBin=ibin;
	  }
      }

    return hist->GetBinCenter(theBin);
  }

bool sortPt (TLorentzVector*a,TLorentzVector*b) { return (a->Pt()>b->Pt()); }

/*
double DeltaR(double eta1, double phi1, double eta2, double phi2) {

double r = sqrt(pow((eta2-eta1),2)+pow((acos(cos(phi2))-acos((phi1))),2));
return r;

}
*/

void printEvent(int i, int nev, time_t t, int total)
{
  if( i%nev==0)
    {
      time_t now=time(NULL);
      float dt= difftime(now,t);
      cout<<setw(10)<<i<<"  processed events"
	  <<"   in "<<setw(6)<<dt<< "s  ==> Rate : "<<i/dt<<" Hz  "<<(i*100)/(total) <<"% "<<endl;
    }
}


// to rerun the official susy analysis
MyAna::MyAna(int dataormc,	
		   TString geometry,
		   float smearing,
		   int nmax,
		   float drcutclus,
		   float deltaTCut, float drcut, int nhits, int cellsize): 
  TupleObjectsFiller(), 
  m_jetCuts(), 
  m_drcutclus(drcutclus),
  m_drcut(drcut),
  m_nmax(nmax),
  m_smearing(smearing),
  m_deltaTCut(deltaTCut),
  m_general("Selection",30),
  m_geometry(geometry),
  m_dataormc(dataormc),
  m_nhits(nhits),
  m_cellsize(cellsize)

{


  //============================================
  // Energy
  //============================================ 

  h_hits_energy = new TH1F("hits_energy","hits_energy",100,0,0.2);
  h_hits_energy->GetXaxis()->SetTitle("energy [MeV]");
  h_hits_energy->GetYaxis()->SetTitle("Nb of hits");


  h_deltaR=new TH1F("deltaR","deltaR",100,0,10);
  h_counterJet1_vs_eta= new TH2F("counterJet1_vs_eta","",17,2.5,4.2,5000,0,5000);
  h_counterJet2_vs_eta = new TH2F("counterJet2_vs_eta","",17,2.5,4.2,5000,0,5000);
  h_counterJet3_vs_eta = new TH2F("counterJet3_vs_eta","",17,2.5,4.2,5000,0,5000);

  h_cell_radius = new TH1F("cell_radius", "cell_radius", 20, 100, 600);
  h_nb_jets = new TH1F("nb_jets", "nb_jets", 50, 0 ,70);
  h_nb_hsjets = new TH1F("nb_hsjets", "nb_hsjets", 20, 0, 20);


  //h_cell_matrix = new TH2F("cell_matrix", "cell_matrix",120,-600, 600, 120, -600, 600);
  
  

      h_nb_cells_perLayer.resize(DetConst::NbLayers+1);//UGLY
      h_nb_cells_en_perLayer.resize(DetConst::NbLayers+1);//UGLY
      h_nb_cells_2hits_perLayer.resize(DetConst::NbLayers+1);//UGLY
      h_nb_cells_dt_perLayer.resize(DetConst::NbLayers+1);
      pr_nb_cells_en_perLayer.resize(DetConst::NbLayers+1);
      pr_nb_cells_2hits_perLayer.resize(DetConst::NbLayers+1);
      pr_nb_cells_dt_perLayer.resize(DetConst::NbLayers+1);
      h_cell_matrix_perLayer.resize(DetConst::NbLayers+1);
      pr_nb_big_cells_en_perLayer.resize(DetConst::NbLayers+1);
      pr_maxnb_big_cells_en_perLayer.resize(DetConst::NbLayers+1);

      for (int i=0;i<h_nb_cells_perLayer.size();i++)
      {

        char tmp[100], tmp1[100], tmp2[100], tmp3[100], tmp4[100], tmp5[100], tmp6[100], tmp7[100], tmp8[100], tmp9[100];
        if((i-4)>=0)
        {
        sprintf(tmp,"%s%d","Nb_cells_Layer_",i-4);
        sprintf(tmp1,"%s%d","Nb_cells_en_Layer_",i-4);
        sprintf(tmp2,"%s%d","Nb_cells_2hits_Layer_",i-4);
        sprintf(tmp3,"%s%d", "Nb_cells_dt_Layer_", i-4);
        sprintf(tmp4,"%s%d", "Pr_cells_en_Layer_", i-4);
        sprintf(tmp5,"%s%d", "Pr_cells_2hits_Layer_", i-4);
        sprintf(tmp6,"%s%d", "Pr_cells_dt_Layer_", i-4);
        sprintf(tmp7,"%s%d", "Cell_Matrix_Layer_", i-4);
        sprintf(tmp8,"%s%d", "Pr_bigcells_en_Layer_", i-4);
        sprintf(tmp9,"%s%d", "Pr_max_bigcells_en_Layer_", i-4);

        }
        else if((i-4)<0)
        {
        sprintf(tmp,"%s%d","Nb_cells_Layer__",abs(i-4));
        sprintf(tmp1,"%s%d","Nb_cells_en_Layer__",abs(i-4));
        sprintf(tmp2,"%s%d","Nb_cells_2hits_Layer__",abs(i-4));
        sprintf(tmp3,"%s%d","Nb_cells_dt_Layer__",abs(i-4));
        sprintf(tmp4,"%s%d", "Pr_cells_en_Layer__", abs(i-4));
        sprintf(tmp5,"%s%d", "Pr_cells_2hits_Layer__", abs(i-4));
        sprintf(tmp6,"%s%d", "Pr_cells_dt_Layer__", abs(i-4));
        sprintf(tmp7,"%s%d", "Cell_Matrix_Layer__", abs(i-4));
        sprintf(tmp8,"%s%d", "Pr_bigcells_en_Layer__", abs(i-4));
        sprintf(tmp9,"%s%d", "Pr_max_bigcells_en_Layer__", abs(i-4));

        }


        h_nb_cells_perLayer[i] = new TH1D(tmp, tmp, 20, 120, 600);
        h_nb_cells_en_perLayer[i]= new TH1D(tmp1,tmp1,20,120,600);
        h_nb_cells_2hits_perLayer[i] = new TH1D(tmp2, tmp2, 20, 120, 600);
        h_nb_cells_dt_perLayer[i] = new TH1D(tmp3, tmp3, 20, 120, 600);
        

        pr_nb_cells_en_perLayer[i] = new TProfile(tmp4, tmp4, 20, 120, 600, 0, 1000);
        pr_nb_cells_2hits_perLayer[i] = new TProfile(tmp5, tmp5, 20, 120, 600, 0, 1000);
        pr_nb_cells_dt_perLayer[i] = new TProfile(tmp6, tmp6, 20, 120, 600, 0, 1000);
        pr_nb_big_cells_en_perLayer[i] = new TProfile(tmp8, tmp8, 48, 120, 600, 0 , 500);
        pr_maxnb_big_cells_en_perLayer[i] = new TProfile(tmp9, tmp9, 48, 120, 600, 0, 500);
        h_cell_matrix_perLayer[i] = new TH2F(tmp7, tmp7,120,-600, 600, 120, -600, 600);
  }

  if(EVTDISPLAY)
    {
      h_evDisplayEtaPhi.resize(DetConst::NbLayers+1);//UGLY 
      h_evDisplay.resize(DetConst::NbLayers+1);//UGLY 

      for (int i=0;i<h_evDisplay.size();i++)
	{
	  char tmp[100];
	  sprintf(tmp,"%s%d","Layer_",i-4);

	  h_evDisplay[i]= new TH2F(tmp,tmp,DetConst::NbCellsInXY,-DetConst::RMax,DetConst::RMax,DetConst::NbCellsInXY,-DetConst::RMax,DetConst::RMax);
	  h_evDisplay[i]->GetXaxis()->SetTitle("x [mm]");
	  h_evDisplay[i]->GetYaxis()->SetTitle("y [mm]");
	  
//	  char tmp2[100];
//	  sprintf(tmp2,"%s%d","Layer_",i-4);
//	  h_evDisplayEtaPhi[i]=DetectorGeometry::createTH2Poly(3500,tmp);
//      cout << "hi " << endl;
	  
	  
	}
    }  

  //============================================
  // 
  //============================================ 
  m_starttime = time(NULL);
  m_nb_events_total=0;


}

MyAna::~MyAna()
{
}
const int cell_width_big = 10;
const int xmin = -600, xmax = 600, ymin = -600, ymax = 600;
const int n_bins_big = (xmax - xmin)/cell_width_big;
//double max_nb_cells[120][120][9]={0};
double max_nb_cells[n_bins_big][n_bins_big][9]={0};
double max_r[n_bins_big][n_bins_big]={0};
double max_z=0;
Bool_t MyAna::Process(Long64_t entry)
{   
const int NbLayers=8;
const int cell_width=m_cellsize;

const int n_bins = (xmax - xmin)/cell_width ;
//double max_nb_cells=0;
//double max_r=0;
//double max_z=0;
//cout << n_bins << endl; 
  //============================================
  // 
  //============================================ 

  m_cellMap.clear();  
  m_cellMap2.clear();   
  //============================================
  // 
  //============================================ 

  gStyle->SetOptStat(0);

  bool doPrint=false;

  if(m_nmax>0 && m_nb_events_total>=m_nmax) return true;

  m_nb_events_total++;
  printEvent(m_nb_events_total,100,m_starttime,this->fChain->GetEntries());


  //============================================
  //init
  //============================================
  float weight = 1.;

  using std::cout;
  using std::endl;
  using std::vector;
  
  //============================================
  // read branches
  //============================================

  GetEntryForSelectedBranches(entry);
  fillall();

  //set seed
  gRandom->SetSeed(eventNumber);

  int incr = 0;
  m_general.increment(weight,incr++,"Begining");

 
  m_var_RunNumber=runNumber;
  m_var_EventNumber=eventNumber;
  m_var_weight=1;

  //============================================
  // vertex
  //============================================

  double time_shift=0;
  //  if(TIMESHIFT)time_shift=gRandom->Gaus(0,TIMESMEARING);
  //  cout<<time_shift<<endl;

  m_Vertex_x=0;
  m_Vertex_y=0;
  m_Vertex_z=0;
  m_Vertex_time=0;
  



//=============================================
// JET CATEGORIZATION
//=============================================

/*
if(JET_CATEGORIZATION)
{
  for(int ijet=0;ijet<AntiKt4EMTopoJets_.size();ijet++)
  {

    if(AntiKt4EMTopoJets_[ijet]->Pt()/1000<30)continue;
    if(AntiKt4EMTopoJets_[ijet]->Pt()/1000>40)continue;
    if(abs(AntiKt4EMTopoJets_[ijet]->Eta())<2.5)continue;
    if(abs(AntiKt4EMTopoJets_[ijet]->Eta())>4.2)continue;
    for(int itruthpart=0;itruthpart<TruthPart_pt.size(); itruthpart++)
    {
        if(TruthPart[itruthpart]->status()!=1) continue;
        if(abs(TruthPart[itruthpart]->pdgId())==13) continue; //remove muon
        if(abs(TruthPart[itruthpart]->pdgId())==12) continue; //remove neutrino
        if(abs(TruthPart[itruthpart]->pdgId())==14) continue; //remove neutrino
        if(abs(TruthPart[itruthpart]->pdgId())==16) continue; //remove neutrino
        double drmin=999999.;
        double dr=AntiKt4EMTopoJets_[itruthjet]->DeltaR(*TruthPart_[itruthpart]);
        if(dr<drmin) drmin=dr;
    }
       //       cout<<drmin<<endl;
       if(drmin<0.4)
    
  }
}
*/







 vector<Jet*> jets;
 vector<Jet*> hsjets;
 for(int ijet=0;ijet<AntiKt4EMTopoJets_.size();ijet++)
    {
      if(AntiKt4EMTopoJets_[ijet]->Pt()/1000<30)continue;
      if(abs(AntiKt4EMTopoJets_[ijet]->Eta())<2.5)continue;
      if(abs(AntiKt4EMTopoJets_[ijet]->Eta())>4.2)continue;

      jets.push_back(AntiKt4EMTopoJets_[ijet]);
      
      double drmin=999999.;
       for(int itruthjet=0;itruthjet<AntiKt4TruthJets_.size();itruthjet++)
	 {
	   double dr=AntiKt4TruthJets_[itruthjet]->DeltaR(*AntiKt4EMTopoJets_[ijet]);

	   if(dr<drmin) drmin=dr;
	 }
       //       cout<<drmin<<endl;
       if(drmin<0.4)       hsjets.push_back(AntiKt4EMTopoJets_[ijet]);

      //      h_jet_pt->Fill(AntiKt4EMTopoJets_[ijet]->Pt()/1000);
    }
    h_nb_jets->Fill(jets.size());
    h_nb_hsjets->Fill(hsjets.size());


  //============================================
  // Loop hits and make cells
  //============================================
  vector<Hit*> allHits=HGTDHits_;
  //h_nbHits->Fill(allHits.size()); 
  for (int ihit=0; ihit<allHits.size();ihit++)
    {  

      double x =DetectorGeometry::sign(allHits[ihit]->x())*( ( fabs(allHits[ihit]->x()) - 1. ) + 0.5 ) * 0.5 ;
      double y= DetectorGeometry::sign(allHits[ihit]->y())*( ( fabs(allHits[ihit]->y()) - 1. ) + 0.5 ) * 0.5 ;
      double z=DetectorGeometry::zCenter((allHits[ihit]->z()+1)*allHits[ihit]->barrel_ec(),m_geometry ); 
      if (z<0)x=-x;
      double r = sqrt(pow(x,2)+pow(y,2));
      double distance0=sqrt(pow(x,2)+
			    pow(y,2)+
			    pow(z,2));

      
      double time=allHits[ihit]->time()+TIMEOFFSET;

      if(fabs(time-distance0/CELERITY)>TIMECUT) continue;

      if (r>=DetConst::RMax) continue;
      if (r<DetConst::RMin) continue;
      
      CellAddress ca;
      DetectorGeometry::IsInsideVar(x,
				 y,
				  z,m_cellsize,
				 ca,m_geometry);

      std::map<CellAddress,Cell>::iterator where = m_cellMap.find(ca);
      if ( where != m_cellMap.end() )
      	{
	  where->second.update(allHits[ihit]->energy(),time);
      	}
      else
      	{
	  Cell cell=Cell(ca,m_cellMap.size()+1);
	  cell.update(allHits[ihit]->energy(),time);

    m_cellMap[ca]=cell;	

    
      	}
    	    

      // h_hits_time_vs_energy->Fill(allHits[ihit]->energy(),time-distance0/CELERITY);
      // h_hits_time->Fill(time-distance0/CELERITY);
      // h_hits_energy->Fill(allHits[ihit]->energy());
      // h_hits_bidim->Fill(allHits[ihit]->x(),
      // 			 allHits[ihit]->y(),
      // 			 allHits[ihit]->energy());
    }//for

  
    if(MB)

    {
  //======================================
  // loop cells to make big 1x1cm2 cells  
  //======================================


  for (std::map<CellAddress,Cell>::iterator it=m_cellMap.begin(); it!=m_cellMap.end(); ++it)
  {
      double energy=it->second.e();
      if(energy<ENERGYCUT) continue;
      double dt=it->second.deltaTMax();
      double xcenter=DetectorGeometry::xCenterVar(it->second.address(), m_cellsize);
      double ycenter=DetectorGeometry::yCenterVar(it->second.address(), m_cellsize);
      double zcenter=DetectorGeometry::zCenter(it->second.address(),m_geometry);
      double phicenter=DetectorGeometry::phiVar(it->second.address(), m_cellsize);
      double etacenter=DetectorGeometry::etaVar(it->second.address(),m_cellsize, m_geometry);
      double r = DetectorGeometry::rCenterVar(it->second.address(), m_cellsize);

      CellAddress cb;
      DetectorGeometry::IsInsideVar(xcenter,
         ycenter,
          zcenter,cell_width_big,
         cb,m_geometry);

      std::map<CellAddress,Cell>::iterator where2 = m_cellMap2.find(cb);
      if ( where2 != m_cellMap2.end() )
        {
    where2->second.update(energy,dt);
        }
      else
        {
    Cell cell2=Cell(cb,m_cellMap2.size()+1);
    cell2.update(energy,dt);

    m_cellMap2[cb]=cell2; 
  }
}







for (std::map<CellAddress,Cell>::iterator it=m_cellMap2.begin(); it!=m_cellMap2.end(); ++it)
{
      double energy=it->second.e();
      double dt=it->second.deltaTMax();
      double nbCells = it->second.nbHits();
     
      double xcenter=DetectorGeometry::xCenterVar(it->second.address(), cell_width_big);
      double ycenter=DetectorGeometry::yCenterVar(it->second.address(), cell_width_big);
      double zcenter=DetectorGeometry::zCenter(it->second.address(),m_geometry);
      double phicenter=DetectorGeometry::phiVar(it->second.address(), cell_width_big);
      double etacenter=DetectorGeometry::etaVar(it->second.address(),cell_width_big, m_geometry);
      double r = DetectorGeometry::rCenterVar(it->second.address(), cell_width_big);
      h_cell_matrix_perLayer[it->second.address().layer()+4]->Fill(xcenter,ycenter, nbCells);
     


      pr_nb_big_cells_en_perLayer[it->second.address().layer()+4]->Fill(r, nbCells);
     

     // if(nbCells>max_nb_cells) {
     //  max_nb_cells = nbCells;
     // pr_maxnb_big_cells_en_perLayer[it->second.address().layer()+4]->
      //  max_r = r;
      //  max_z = it->second.address().layer()+4;
      //}
}

}

 double cell_x_big[n_bins_big]={0}, cell_y_big[n_bins_big]={0};
 for(int i=0; i<n_bins_big; i++)
 {
  cell_x_big[i]=(-600+cell_width_big*i)+cell_width_big/2.;
  cell_y_big[i]=(-600+cell_width_big*i)+cell_width_big/2.;
 }


for(int l=0; l<h_nb_cells_perLayer.size(); l++)
{
  for(int i=0; i<n_bins_big; i++)
  {
    for(int j=0; j<n_bins_big; j++)
    {
      if((h_cell_matrix_perLayer[l]->GetBinContent(i+1, j+1))>max_nb_cells[i][j][l])
      {
        max_nb_cells[i][j][l]=h_cell_matrix_perLayer[l]->GetBinContent(i+1, j+1);
        //cout << cell_x_big[i] << "   " << cell_y_big[] << endl;
        max_r[i][j] = sqrt(pow(cell_x_big[i],2)+pow(cell_y_big[j],2));
      }
      
      //if(h_cell_matrix_perLayer[l]->GetBinContent(i+1, j+1)!=0)
    }
  }
}

//for(int i=0; i<n_bins_big; i++)
//{
//  for(int j=0; j<n_bins_big; j++)
//  {
//    cout << max_nb_cells[i][j] << "   " << max_r << endl;
//  }
//}


//pr_maxnb_big_cells_en_perLayer[max_z]->Fill(max_r, max_nb_cells);  
  //============================================
  //Fill Cell Matrix
  //============================================

 double cell_x[n_bins]={0}, cell_y[n_bins]={0};
 for(int i=0; i<n_bins; i++) 
 {
 cell_x[i]=(-600+cell_width*i)+cell_width/2.;
 cell_y[i]=(-600+cell_width*i)+cell_width/2.;
 //cout << cell_x[i] << "  " << cell_y[i] << endl;
 }
 

 //int n_bins_big = (xmax - xmin)/cell_width_big ;


//cout << cell_x_big[0] << "   " << cell_y_big[0];
  //============================================
  // Cell container analysis
  //============================================

  int counter=0;
  int counter1=0;
  int counter2=0;
  int counterHits=0;
  for (std::map<CellAddress,Cell>::iterator it=m_cellMap.begin(); it!=m_cellMap.end(); ++it)
    {
      double energy=it->second.e();
      double dt=it->second.deltaTMax();
      
      if(it->second.nbHits()==0) counter1++;
      
      counter++;
      counterHits+=it->second.nbHits();
      double xcenter=DetectorGeometry::xCenterVar(it->second.address(), m_cellsize);
      double ycenter=DetectorGeometry::yCenterVar(it->second.address(), m_cellsize);
      double zcenter=DetectorGeometry::zCenter(it->second.address(),m_geometry);
      double phicenter=DetectorGeometry::phiVar(it->second.address(), m_cellsize);
      double etacenter=DetectorGeometry::etaVar(it->second.address(),m_cellsize, m_geometry);
      double r = DetectorGeometry::rCenterVar(it->second.address(), m_cellsize);
	    if(energy>ENERGYCUT)
      { 
        
        h_nb_cells_en_perLayer[it->second.address().layer()+4]->Fill(r);
        if(it->second.nbHits()>=2) 
        {
          h_nb_cells_2hits_perLayer[it->second.address().layer()+4]->Fill(r);
          counter2++;
          if(dt>=0.02)
          {
            h_nb_cells_dt_perLayer[it->second.address().layer()+4]->Fill(r);
          }
        }
      //cout << counter2 << "    " << r << endl;
      }
      


      if(EVTDISPLAY)
	{
	  h_evDisplay[it->second.address().layer()+4]->Fill(xcenter,ycenter,energy);
	  //h_evDisplayEtaPhi[it->second.address().layer()+4]->Fill(abs(etacenter),phicenter,energy);
	}
	  
      double distance=sqrt(pow(xcenter-m_Vertex_x,2)+
			   pow(ycenter-m_Vertex_y,2)+
			   pow(zcenter-m_Vertex_z,2));
	  
      double time=it->second.time(m_smearing)-distance/CELERITY;//ns
      

	  //	  h_cell_energy->Fill(energy);  
	  //	  h_cell_energy_perLayer[abs(it->second.address().layer())-1]->Fill(energy);  
	
	}
  //  cout<<counter<<"// "<<counterHits/float(counter)<<"  "<<counter1<<" "<<counter2/float(counter1)<<endl;



      
  

  for(int l=0; l<h_nb_cells_perLayer.size(); l++) 
  {
  for(int i=0; i<n_bins; i++)
  {
  for(int j=0; j<n_bins; j++)
  {
    double radius = sqrt(pow(cell_x[i],2)+pow(cell_y[j],2));
    if(radius > 120 && radius < 600)
    h_nb_cells_perLayer[l]->Fill(radius);
 
  }
}
}

  


  //============================================
  // Loop over jets
  //============================================

  //  cout<<"======"<<endl;
 for(int ijet=0;ijet<hsjets.size();ijet++)
    {


      double eta_jet = hsjets.at(ijet)->Eta();
      double phi_jet = hsjets.at(ijet)->Phi();
      double theta_jet = 2*atan(exp(-eta_jet));
      int counterJet1=0;
      int counterJet2=0;
      int counterJet3=0;
      double counter_cells[DetConst::NbLayers+1] = {0}, counter_hits[DetConst::NbLayers+1] = {0}, counter_time[DetConst::NbLayers+1] ={0};
      for (std::map<CellAddress,Cell>::iterator it=m_cellMap.begin(); it!=m_cellMap.end(); ++it)
	    {
        

        
	      double deltaR = DetectorGeometry::deltaRVar(it->second.address(),hsjets[ijet]->Eta(),hsjets[ijet]->Phi(), m_cellsize,m_geometry);
	      //cout << it->second.address().layer() << "   " << deltaR <<endl;
        h_deltaR->Fill(deltaR);
          if(it->second.e()<ENERGYCUT) continue;
    
              if(deltaR<m_drcut) 
              {
                counter_cells[it->second.address().layer()+4]++;
                //double r_jet = fabs(DetectorGeometry::zCenter(it->second.address(), m_geometry)*sin(theta_jet)/cos(theta_jet));
                //cout << DetectorGeometry::zCenter(it->second.address(), m_geometry) << endl;
                if(it->second.nbHits()>=m_nhits)
                {
                  counter_hits[it->second.address().layer()+4]++;
                  if(it->second.deltaTMax()>=0.02)
                  {
                    counter_time[it->second.address().layer()+4]++;
                  }
                }
              }
	      
	  //cout<<it->second.nbHits()<<"  "<<it->second.deltaTMax()<<endl;

	    }
      for(int i=0; i<DetConst::NbLayers+1; i++)
      {
        //cout << counter_cells[i] << endl;
        //cout << DetectorGeometry::zCenter(i-4, m_geometry) << endl;
      double r_jet = fabs(DetectorGeometry::zCenter(i-4, m_geometry)*sin(theta_jet)/cos(theta_jet));
      pr_nb_cells_en_perLayer[i]->SetBit(TH1::kIsAverage);
      pr_nb_cells_en_perLayer[i]->Fill(r_jet, counter_cells[i]);

      pr_nb_cells_2hits_perLayer[i]->SetBit(TH1::kIsAverage);
      pr_nb_cells_2hits_perLayer[i]->Fill(r_jet, counter_hits[i]);
      
      pr_nb_cells_dt_perLayer[i]->SetBit(TH1::kIsAverage);
      pr_nb_cells_dt_perLayer[i]->Fill(r_jet, counter_time[i]);
      h_counterJet1_vs_eta->Fill(abs(hsjets[ijet]->Eta()),counterJet1);
      h_counterJet2_vs_eta->Fill(abs(hsjets[ijet]->Eta()),counterJet2);
      h_counterJet3_vs_eta->Fill(abs(hsjets[ijet]->Eta()),counterJet3);
    }
    }





for(int i=0; i<h_evDisplay.size(); i++) 
{
  //h_evDisplay.at(i)->Reset();
}

jets.clear();
hsjets.clear();

for(int i=0; i<h_cell_matrix_perLayer.size(); i++)
{
  h_cell_matrix_perLayer[i]->Reset();

}

return kTRUE       ;
      
}




void MyAna::Terminate()
{




for(int l=0; l<9; l++)
{
  for(int i=0; i<n_bins_big; i++)
  {
    for(int j=0; j<n_bins_big; j++)
    {
      
pr_maxnb_big_cells_en_perLayer[l]->Fill(max_r[i][j], max_nb_cells[i][j][l]); 
      //if(max_nb_cells[i][j][l]>25)
      //cout << max_nb_cells[i][j][l] << "   " << max_r[i][j] << endl;
      

    }
  }
}
  

  

  int di=8;
  for(int i=0; i<4; i++)
  { 
    h_nb_cells_perLayer.at(i)->Add(h_nb_cells_perLayer.at(i+di));
    h_nb_cells_en_perLayer.at(i)->Add(h_nb_cells_en_perLayer.at(i+di));   
    h_nb_cells_2hits_perLayer.at(i)->Add(h_nb_cells_2hits_perLayer.at(i+di));
    h_nb_cells_dt_perLayer.at(i)->Add(h_nb_cells_dt_perLayer.at(i+di));


    pr_nb_cells_en_perLayer.at(i)->Add(pr_nb_cells_en_perLayer.at(i+di));   
    pr_nb_cells_2hits_perLayer.at(i)->Add(pr_nb_cells_2hits_perLayer.at(i+di));
    pr_nb_cells_dt_perLayer.at(i)->Add(pr_nb_cells_dt_perLayer.at(i+di));
    
    pr_nb_big_cells_en_perLayer.at(i)->Add(pr_nb_big_cells_en_perLayer.at(i+di));
    pr_maxnb_big_cells_en_perLayer.at(i)->Add(pr_maxnb_big_cells_en_perLayer.at(i+di));

    h_nb_cells_2hits_perLayer.at(i)->Divide( h_nb_cells_perLayer.at(i) );

    h_nb_cells_en_perLayer.at(i)->Divide( h_nb_cells_perLayer.at(i) );
    h_nb_cells_dt_perLayer.at(i)->Divide( h_nb_cells_perLayer.at(i) );

    h_nb_cells_en_perLayer.at(i)->Scale(100); 
    h_nb_cells_2hits_perLayer.at(i)->Scale(100); 
    h_nb_cells_dt_perLayer.at(i)->Scale(100); 

    h_nb_cells_en_perLayer.at(i)->SetMinimum(0.2);
    h_nb_cells_en_perLayer.at(i)->SetMaximum(100);

    h_nb_cells_2hits_perLayer.at(i)->SetMinimum(0.02);
    h_nb_cells_2hits_perLayer.at(i)->SetMaximum(20);

    h_nb_cells_dt_perLayer.at(i)->SetMinimum(0.002);
    h_nb_cells_dt_perLayer.at(i)->SetMaximum(10);

    di=di-2;
  }


  using std::cout;
  using std::endl;
  cout<<m_general<<endl;

}


 

int MyAna::getRBins(double r)
{
  int irbin=-1;

  for (int ir=0;ir<m_rbins.size();ir++)
    {
      //      cout<<m_rbins[ir]<<"  "<<m_rbins[ir+1]<<"   "<< r<<"  "<<bool(m_rbins[ir]>r)<<"  "<<bool(r<m_rbins[ir+1])  <<endl;
      if(m_rbins[ir]<r && r<m_rbins[ir+1] )
	{
	  irbin=ir;	  
	}
    }
  //  cout<<r<<" _ "<<irbin<<endl;
  return irbin;
}






























//======================================================
//======================================================
//======================================================
//toto



void MyAna::doCellLevelJetAnalysis(int ijet)
{

  MyJet* jet= m_jets[ijet];

  //=============================================
  // Make time plot
  //=============================================
  h_time_tmp->Reset();

  std::vector <Cell> cellVec=jet->cellVec();
  for (int iCell=0;iCell<cellVec.size();iCell++)
    {
      Cell cell = cellVec[iCell];
      double distance=DetectorGeometry::distance(cell.address(),m_Vertex_x,m_Vertex_y,m_Vertex_z,m_geometry);
      double time=cell.time(m_smearing)-distance/CELERITY;
      double clusEnergy=m_clusters[cell.clusIndex()]->E();
      h_time_tmp->Fill(time,clusEnergy);

    }//loop cells
  


  //=============================================
  // compute time var
  //=============================================

  float timeAtMax= slidingWindow(h_time_tmp,WIDTHWINDOW);//width should be even

  float ITF_den=0;
  float ITF_num=0;
  float eta=0;
  float phi=0;
  float sumE=0;
  for (int iCell=0;iCell<cellVec.size();iCell++)
    {
      Cell cell = cellVec[iCell];
      double distance=DetectorGeometry::distance(cell.address(),m_Vertex_x,m_Vertex_y,m_Vertex_z,m_geometry);
      double time=cell.time(m_smearing)-distance/CELERITY;
      double clusEnergy=m_clusters[cell.clusIndex()]->E();
      double weight=clusEnergy;
      double deltaT=time-timeAtMax;
      
      if(jet->isHS())
	{
	  h_jet_HGTDdeltaT_hs->Fill(deltaT,weight);
	}
      else 
	{
	  h_jet_HGTDdeltaT_pu->Fill(deltaT,weight);
	}
      

      ITF_den+=weight;//clusEnergy;
      if(fabs(deltaT)<m_deltaTCut)
	{
	  sumE+=clusEnergy;
	  ITF_num+=weight;//clusEnergy;
	  double etaCell=DetectorGeometry::eta(cell.address(),m_geometry);
	  double phiCell=DetectorGeometry::phi(cell.address());
	  eta+=etaCell*clusEnergy;
	  phi+=phiCell*clusEnergy;
	}
    }//loop cells

  float ITF=1;
  if (ITF_den!=0) ITF=ITF_num/ITF_den;
  jet->setITF(ITF);
  if(sumE!=0)
    {
      eta=eta/sumE;
      phi=phi/sumE;
    }
	  



  //=============================================
  // compute width
  //=============================================
    float   width=0;
  for (int iCell=0;iCell<cellVec.size();iCell++)
    {
      Cell cell = cellVec[iCell];
      double distance=DetectorGeometry::distance(cell.address(),m_Vertex_x,m_Vertex_y,m_Vertex_z,m_geometry);
      double time=cell.time(m_smearing)-distance/CELERITY;
      double clusEnergy=m_clusters[cell.clusIndex()]->E();
      double deltaT=time-timeAtMax;
      if(fabs(deltaT)<m_deltaTCut)
	{

	  double dr= DetectorGeometry::deltaR(cell.address(),eta,phi,m_geometry);
	  width+=dr*clusEnergy;
	}
    }
  if(sumE!=0) width=width/sumE;
  jet->setWidth(width);


  //=============================================
  // draw time
  //=============================================
  if(MAKETIMEPLOT)
    {
      TCanvas* canvasTMP=0;
      canvasTMP=new TCanvas("canvasTMP","canvasTMP",1000,600);
      
      h_time_tmp->Draw();
      
      TLatex latex=TLatex();
      latex.SetTextSize(0.02);
      latex.SetNDC();
	  
	  
      TLine line = TLine(m_TrueVertex_time/CELERITY,0,m_TrueVertex_time/CELERITY,h_time_tmp->GetMaximum()*2);
      line.SetLineColor(2);
      line.Draw();
	  
      TLine line2 = TLine(timeAtMax,0,timeAtMax,h_time_tmp->GetMaximum()*2);
      line2.SetLineColor(4);
      line2.Draw();

      char tmp[100];
      sprintf(tmp,"%s%d%s%d%s","TIME_",m_nb_events_total,"_",ijet,"bis.gif");
      if(MAKETIMEPLOT)canvasTMP->Print(tmp);
	      

    }



}
