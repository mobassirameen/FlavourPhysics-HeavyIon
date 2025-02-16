{
//=========Macro generated from canvas: c1_n2/A canvas
//=========  (Sat Jan  6 16:08:14 2018) by ROOT version5.34/30
   TCanvas *c1_n2 = new TCanvas("c1_n2", "A canvas",201,117,1050,963);
   c1_n2->Range(0,0,140,60);
   c1_n2->SetFillColor(0);
   c1_n2->SetBorderMode(0);
   c1_n2->SetBorderSize(2);
   c1_n2->SetFrameBorderMode(0);
   TArrow *arrow = new TArrow(36,42,94,42,0.02,"->-");
   arrow->SetFillColor(1);
   arrow->SetFillStyle(1001);
   arrow->SetLineWidth(3);
   arrow->Draw();
   arrow = new TArrow(36,30,94,30,0.02,"-<-");
   arrow->SetFillColor(1);
   arrow->SetFillStyle(1001);
   arrow->SetLineWidth(3);
   arrow->Draw();
   TCurlyLine *curlyline = new TCurlyLine(62.10325,30.07389,80,20,0.02,0.01);
   curlyline->SetLineWidth(3);
   curlyline->Draw();
   arrow = new TArrow(80,20,94,25,0.02,"-<-");
   arrow->SetFillColor(1);
   arrow->SetFillStyle(1001);
   arrow->SetLineWidth(3);
   arrow->Draw();
   
   TArc *arc = new TArc(98.10707,43.44828,0,0,360);
   arc->Draw();
   
   TEllipse *ellipse = new TEllipse(95.96558,29.40887,6.022945,15.14778,0,360,0);
   ellipse->SetLineWidth(3);
   ellipse->Draw();
   
   ellipse = new TEllipse(33.32696,35.39409,4.550669,9.08867,0,360,0);
   ellipse->SetLineWidth(3);
   ellipse->Draw();
   arrow = new TArrow(80,20,94,15,0.02,"->-");
   arrow->SetFillColor(1);
   arrow->SetFillStyle(1001);
   arrow->SetLineWidth(3);
   arrow->Draw();
   TCurlyArc *curlyarc = new TCurlyArc(62.12851,30,10.92997,180,180,0.02,0.01);
   curlyarc->Draw();
   curlyarc = new TCurlyArc(62,30,11,0,180,0.013,0.013);
   curlyarc->SetWavy();
   curlyarc->SetLineWidth(4);
   curlyarc->Draw();
   TLatex *   tex = new TLatex(95.58233,28.69266,"f");
   tex->SetTextFont(12);
   tex->SetTextSize(0.05733945);
   tex->SetTextAngle(3.576334);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(38.65462,42.31651,"d or s");
   tex->SetTextFont(42);
   tex->SetTextSize(0.03440367);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(30.36145,38.32569,"B^{0}");
   tex->SetTextFont(42);
   tex->SetTextSize(0.04357798);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(31.06426,34.6789,"or");
   tex->SetTextFont(42);
   tex->SetTextSize(0.04357798);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(30.36145,30.06881,"B_{s}");
   tex->SetTextFont(42);
   tex->SetTextSize(0.04357798);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(41.3253,30.2,"#bar{b}");
   tex->SetTextFont(42);
   tex->SetTextSize(0.03440367);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(78.71486,30.8945,"#bar{q}'");
   tex->SetTextFont(42);
   tex->SetTextSize(0.03440367);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(58.19277,30.66,"#bar{q}^{u}");
   tex->SetTextFont(42);
   tex->SetTextSize(0.03440367);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(45.54217,26.2156,"V_{q}^{*}u_{b}");
   tex->SetTextFont(42);
   tex->SetTextSize(0.04357798);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(70.28112,27.93578,"V_{q}u_{q'}");
   tex->SetTextFont(42);
   tex->SetTextSize(0.04357798);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(83.58209,23.08789,"#bar{q}");
   tex->SetTextFont(42);
   tex->SetTextSize(0.03440367);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(82.77842,16.60333,"q");
   tex->SetTextFont(42);
   tex->SetTextSize(0.03440367);
   tex->SetLineWidth(2);
   tex->Draw();
   c1_n2->Modified();
   c1_n2->cd();
   c1_n2->SetSelected(c1_n2);
   c1_n2->ToggleToolBar();
}
