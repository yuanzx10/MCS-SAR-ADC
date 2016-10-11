#include<iostream>
#include<cmath>
#include<vector>
#include"TMath.h"
using namespace std;

void set_Graph(TGraph* gr,const char* xtitle,const char* ytitle);
double get_binsize(int M, int L, vector<double>& binsize);

void mismatch_effect_mc(){
    int M = 5, L= 4;
    int codes = pow(2.0,M+L+1);
    vector<double> binsize(codes);


    int MC_run = 400;
    vector<vector<double>> DNL(codes, vector<double>(MC_run));
    vector<vector<double>> INL(codes, vector<double>(MC_run));
    TGraph *gr_DNL = new TGraph(codes*MC_run); gr_DNL->SetTitle("DNL Plot");
    TGraph *gr_INL = new TGraph(codes*MC_run); gr_INL->SetTitle("DNL Plot");
    for(int iMC=0;iMC<MC_run;iMC++){
        cout<<iMC<<endl;
        double LSB = get_binsize(M,L,binsize);

        // calculate the DNL, we have force DNL[0]=0,DNL[codes-1]=0
        for(int i=0;i<codes;i++){
            DNL[i][iMC]=binsize[i]/LSB-1;
            gr_DNL->SetPoint(i+iMC*codes,i,DNL[i][iMC]);
        }
        // calculate the INL
        INL[0][iMC] = 0; gr_INL->SetPoint(iMC*codes,0,0);
        for(int i=1;i<codes;i++){
            INL[i][iMC] = DNL[i-1][iMC]+INL[i-1][iMC];
            gr_INL->SetPoint(i+iMC*codes,i,INL[i][iMC]);
        }

    }// the MC loop

    // plot the result of MC simulation
    TCanvas *c1 = new TCanvas("c1","DNL and INL",0,0,800,600);
	c1->Divide(1,2);
	c1->cd(1);
	gr_DNL->Draw("A*");	set_Graph(gr_DNL,"Code","LSB");
	c1->cd(2);
	gr_INL->Draw("A*");	set_Graph(gr_INL,"Code","LSB");

	// calculate the mean and std value
	TGraph *gr_std_DNL = new TGraph(codes);	
	gr_std_DNL->SetTitle("Standard Deviation of DNL");
	TGraph *gr_std_INL = new TGraph(codes);
	gr_std_INL->SetTitle("Standard Deviation of INL");

	double mean_dnl, stdv_dnl, max_dnl=0.0;
	vector<double> mean_DNL(codes);
	vector<double> stdv_DNL(codes);
	double mean_inl, stdv_inl, max_inl=0.0;
	vector<double> mean_INL(codes);
	vector<double> stdv_INL(codes);
	for(int i=0;i<codes;i++){ // get mean
		mean_dnl = 0.; mean_inl = 0.0;
		for(int j=0;j<MC_run;j++){mean_dnl += DNL[i][j];mean_inl += INL[i][j];}
		mean_dnl /= MC_run;	mean_inl /= MC_run;
		mean_DNL[i] = mean_dnl; mean_INL[i] = mean_inl;
	}
	for(int i=0;i<codes;i++){
		stdv_dnl = 0.0; stdv_inl = 0.0;
		for(int j=0;j<MC_run;j++){
			stdv_dnl += (DNL[i][j]-mean_DNL[i])*(DNL[i][j]-mean_DNL[i]);
			stdv_inl += (INL[i][j]-mean_INL[i])*(INL[i][j]-mean_INL[i]);
		}
		stdv_dnl = sqrt(stdv_dnl/MC_run);	stdv_inl = sqrt(stdv_inl/MC_run);
		stdv_DNL[i] = stdv_dnl;	stdv_INL[i] = stdv_inl;
		max_dnl = (stdv_dnl>max_dnl)?stdv_dnl:max_dnl; // get maximum std value
		max_inl = (stdv_inl>max_inl)?stdv_inl:max_inl;
		gr_std_DNL->SetPoint(i,i,stdv_dnl);
		gr_std_INL->SetPoint(i,i,stdv_inl);
	}
	for(int i=0;i<codes;i=i+2){
		cout<<i<<"\t"<<mean_DNL[i]<<"\t"<<stdv_DNL[i]
		cout<<"\t"<<mean_INL[i]<<"\t"<<stdv_INL[i]<<endl;
	}
    cout<<endl<<"maximum DNL std = \t"<<max_dnl;
    cout<<"\t maximum DNL std = \t"<<max_inl<<endl;
    cout<<"peridicted DNL std = \t"<<pow(2.,L)*sqrt(pow(2.,M)-1)*0.005<<endl;


	TCanvas *c2 = new TCanvas("c2","DNL and INL",800,0,800,600);
	c2->Divide(1,2);
	c2->cd(1);
	gr_std_DNL->Draw("A*");set_Graph(gr_std_DNL,"Code","LSB");
	c2->cd(2);
	gr_std_INL->Draw("A*");set_Graph(gr_std_INL,"Code","LSB");


	TCanvas *c3 = new TCanvas("c3","Histogram of DNL",0,550,800,600);
	TH1D *h1 = new TH1D("c1","Maximum DNL Distribution Plot",25,-1,1);
	for(int i=0;i<MC_run;i++){h1->Fill(DNL[codes/2-1][i]);}
	h1->Draw();
	h1->Fit("gaus");
	h1->GetXaxis()->SetTitle("LSB");
	h1->GetYaxis()->SetTitle("Count");

}


void set_Graph(TGraph* gr,const char* xtitle,const char* ytitle){
	// set X axis
	double x_max = gr->GetX()[gr->GetN()-1];
	double x_min = gr->GetX()[0];
	gr->GetXaxis()->SetLimits(x_min,x_max);
	gr->GetXaxis()->SetTitle(xtitle);
	gr->GetXaxis()->CenterTitle();
	gr->GetXaxis()->SetTitleFont(62);
	gr->GetXaxis()->SetTitleSize(0.05);
	gr->GetXaxis()->SetLabelFont(62);
	gr->GetXaxis()->SetLabelSize(0.05);
	// set y axis
	gr->GetYaxis()->SetTitle(ytitle);
	gr->GetYaxis()->CenterTitle();
	gr->GetYaxis()->SetTitleFont(62);
	gr->GetYaxis()->SetTitleSize(0.05);
	gr->GetYaxis()->SetLabelFont(62);
	gr->GetYaxis()->SetLabelSize(0.05);
	//gr->Draw("APL");
}

double get_binsize(int M, int L, vector<double>& binsize){
    vector<double> cap(M+L);    // postive capacitor array
    vector<double> can(M+L);    // negative capacitor array
    vector<int> D(M+L+1);       // Bits
    int codes = pow(2.0,M+L+1); // overall codes

    double vrefp = 1.0, vrefn = 0.0, vcm = 0.5; // reference voltages

    double sigma = 0.005;        // variation of the unit capacitor
    TRandom *r = new TRandom(0);// random number generator
    //m_rand->SetSeed(3);

    // define the capacitors including the parasitic
    double cbp = 1.0 + sigma*r->Gaus(0.0,1.0); // bridge capacitor
    double cbn = 1.0 + sigma*r->Gaus(0.0,1.0);
    double cmp,cmn,clp,cln; // cm, cl is the total capacitance
    for(int i=0;i<M+L;i++){
        if(i<L){ // LSB capacitor
            cap[i] = pow(2.0,i)+sqrt(pow(2.0,i))*sigma*r->Gaus(0.0,1.0);
            clp += cap[i];
            can[i] = pow(2.0,i)+sqrt(pow(2.0,i))*sigma*r->Gaus(0.0,1.0);
            cln += can[i];
        }
        else{   // MSB capacitors
            cap[i] = pow(2.0,i-L)+sqrt(pow(2.0,i-L))*sigma*r->Gaus(0.0,1.0);
            cmp += cap[i];
            can[i] = pow(2.0,i-L)+sqrt(pow(2.0,i-L))*sigma*r->Gaus(0.0,1.0);
            cmn += can[i];
        }
    }
    // add the parasitic capacitance here if you want, add it to cm,cl

    // coef used to calculate the DAC voltage change
    double p_m2m = (cbp+clp)/(cmp*(cbp+clp)+cbp*clp);
    double p_l2m = (cbp+0.0)/(cmp*(cbp+clp)+cbp*clp);
    double n_m2m = (cbn+cln)/(cmn*(cbn+cln)+cbn*cln);
    double n_l2m = (cbn+0.0)/(cmn*(cbn+cln)+cbn*cln);

    // calculate the voltage change due to bit conversion
    vector<double> p_dac(M+L);  //voltage change when postive capacitor switch 
    vector<double> n_dac(M+L);  
    for(int i=0;i<M+L;i++){
        if(i<L){    // LSB capacitor conversion
            p_dac[i] = (vrefp-vcm)*cap[i]*p_l2m + (vcm-vrefn)*can[i]*n_l2m;
            n_dac[i] = (vcm-vrefn)*cap[i]*p_l2m + (vrefp-vcm)*can[i]*n_l2m;
        }
        else{       // MSB capacitor conversion
            p_dac[i] = (vrefp-vcm)*cap[i]*p_m2m + (vcm-vrefn)*can[i]*n_m2m;
            n_dac[i] = (vcm-vrefn)*cap[i]*p_m2m + (vrefp-vcm)*can[i]*n_m2m;
        }
    }

    // for loop fpr all code from 1 to codes-2
    double totalbin = 0;
    for(int icode=1;icode<codes-1;icode++){
        // get the binary of icode
        int num = icode;
        for(int i=0;i<=M+L;i++){
            D[i] = 0;
            if(num%2==1) {D[i]=1;}
            num = num/2;
        }

        // calculate the bin size of icode
        // this is according to the algorithm of MCS SAR
        double dlimit=-50000., ulimit = 50000.;
        double vi = 0.0;
        for(int j=M+L;j>=0;j--){
            if(D[j]>0){
                if(vi>dlimit) {dlimit = vi;}
                if(j>0) {vi = vi + n_dac[j-1];}
            }
            else{
                if(vi<ulimit) {ulimit = vi;}
                if(j>0) {vi = vi - p_dac[j-1];}
            }
        }
        binsize[icode] = ulimit-dlimit;
        if(binsize[icode]<0) {binsize[icode]=0;} //if binsize < 0, missing code
        totalbin += binsize[icode];
    } // end of icode loop
    double LSB = totalbin/(codes-2);
    binsize[0] = LSB; binsize[codes-1] = LSB;

    return LSB;

}
