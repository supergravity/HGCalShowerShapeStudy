void plotRatio(int npoints, TGraphErrors *gd, TGraphErrors *gs, string name, string title){



  
    int W = 800;
    int H = 600;
    
    int H_ref = 600; 
    int W_ref = 800; 
    
    // references for T, B, L, R
    float T = 0.08*H_ref;
    float B = 0.12*H_ref; 
    float L = 0.12*W_ref;
    float R = 0.04*W_ref;
    
    TCanvas* c = new TCanvas("c","c",50,50,W,H);
    
    c->SetLeftMargin( L/W );
    c->SetRightMargin( R/W );
    c->SetTopMargin( T/H );
    c->SetBottomMargin( B/H );
    
    gStyle->SetOptStat(0);
    TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.2,1.0,1.0,21);
    TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.198,22);
    pad1->SetFillColor(0);
    pad2->SetFillColor(0);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();

    gs->GetYaxis()->SetTitle(title.c_str());
    gs->GetXaxis()->SetTitle("Beam Energy");

    gd->GetYaxis()->SetTitle("<Energy> [MIPS]");
    gd->GetXaxis()->SetTitle("Beam Energy");

    
    gs->Draw("AP");
    gd->Draw("Psame");

    //TLegend *leg = new TLegend(0.7669173,0.7582609,0.9461153,0.9147826,NULL,"brNDC");
    TLegend *leg = new TLegend(0.1566416,0.7326087,0.3358396,0.8891304,NULL,"brNDC");
    leg->SetBorderSize(0);
    ///for LINEAR SCALE
    leg->SetTextFont(62);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(1001);
    
    leg->AddEntry(gd,"Data","LP");
    leg->AddEntry(gs,"Simulation","LP");
    leg->Draw();
    
    pad2->cd();
    

    

    
    double ratio[100];
    double x[100];
    
    for(int ii=0; ii<gd->GetN(); ii++){
      
      double stmpy, stmpx, dtmpx, dtmpy;
      gd->GetPoint(ii,dtmpx,dtmpy);
      gs->GetPoint(ii,stmpx,stmpy);
    
      cout<<"i : x : y "<<ii<<" "<<dtmpx<<" "<<dtmpy<<endl;
      ratio[ii] = dtmpy/stmpy;
      x[ii] = dtmpx;

      
    }
    
    TGraph *hratio = new TGraph(gd->GetN(),x,ratio);
    



    hratio->GetXaxis()->SetLabelSize(0.11);
    hratio->GetYaxis()->SetLabelSize(0.11);
    hratio->GetYaxis()->SetTitleSize(0.09);
    
    hratio->GetXaxis()->SetLabelFont(42);
    hratio->GetXaxis()->SetLabelSize(0.11);
    hratio->GetXaxis()->SetTitleSize(0.035);
    hratio->GetXaxis()->SetTitleFont(42);
    hratio->GetYaxis()->SetTitle("#frac{Data}{SM}");
    hratio->GetYaxis()->SetLabelFont(42);
    hratio->GetYaxis()->SetLabelSize(0.11);
    hratio->GetYaxis()->SetTitleSize(0.13);
    hratio->GetYaxis()->SetTitleOffset(0.37);
    hratio->SetMarkerColor(4);
    hratio->SetMarkerStyle(20);

    hratio->GetYaxis()->SetTickLength(0.01);
    
    hratio->GetYaxis()->SetTitle("Data/SM");
    hratio->SetMaximum(1.15);
    //hratio->SetMinimum(0);
    hratio->SetMinimum(0.8);
    hratio->Draw("AP");
    
    
    //TLine *l = new TLine(xlow,1.,xhigh,1.);
    TLine *l = new TLine(x[0],1.,x[gd->GetN()-1],1.);
    l->SetLineColor(2);
    l->SetLineStyle(2);
    l->SetLineWidth(2);
    l->Draw("sames");
    

   c->Modified();
    c->Update();
    
    c->Print(Form("plots/%s.png",name.c_str()));
    c->Print(Form("plots/%s.C",name.c_str()));
    c->Print(Form("plots/%s.pdf",name.c_str()));

 
}
