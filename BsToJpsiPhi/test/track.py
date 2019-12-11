#
# Original Author:Julie Hogan (FNAL)
# edite by: Jhovanny Andres Mejia Guisao
#        

from DataFormats.FWLite import Handle, Events
import ROOT

events = Events('root://cms-xrd-global.cern.ch//store/data/Run2017B/Charmonium/MINIAOD/PromptReco-v1/000/297/046/00000/88DF6C6A-4556-E711-A5C0-02163E01A630.root')
tracks = Handle("std::vector<pat::PackedCandidate>")

histogram = ROOT.TH1F("histogram", "histogram", 100, 0, 8)

i = 0
for event in events:
    print "Event", i
    event.getByLabel("packedPFCandidates", tracks)
    j = 0
    for track in tracks.product():
        print "    Track", j, track.charge() / track.pt(), track.phi(), track.eta(), track.dxy(), track.dz()
        histogram.Fill(track.pt())
        j += 1
    i += 1
    if i >= 5: break


c = ROOT.TCanvas ( "c" , "c" , 800, 800 )
c.cd()
histogram.Draw()
c.SaveAs("PtFW.png")

