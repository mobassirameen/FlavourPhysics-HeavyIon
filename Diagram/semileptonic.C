{
//=========Macro generated from canvas: c1/A canvas
//=========  (Sat Jan  6 13:51:05 2018) by ROOT version5.34/30
   TCanvas *c1 = new TCanvas("c1", "A canvas",745,92,825,899);
   c1->Range(0,0,140,60);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
   TLine *line = new TLine(32,48,116,48);
   line->SetLineWidth(3);
   line->Draw();
   line = new TLine(32,32,116,32);
   line->SetLineWidth(3);
   line->Draw();
   
  
      tex = new TLatex(97,40,"W^{+}");
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(91.77258,49.7561,"V^{*}_{ts}");
   tex->SetLineWidth(2);
   tex->Draw();
   
     
      tex = new TLatex(38,33,"s");
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(100,33,"b");
   tex->SetLineWidth(2);
   tex->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
   c1->ToggleToolBar();
}
