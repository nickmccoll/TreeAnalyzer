
{
  for(float mass = 800; mass<=3500; mass+=50){
    TString line = TString::Format("/Users/nmccoll/Dropbox/Work/GitRepositories/TreeAnalyzer/framework/HiggsAnalysis/CombinedLimit/exe/combine -M Significance  -m %f ../combined.root",mass);
      gSystem->Exec(line.Data());
  }
}

{
  TGraph * g = new TGraph();
     g->SetPoint( 0 ,  800, RooStats::SignificanceToPValue(        0));
     g->SetPoint( 1 ,  850, RooStats::SignificanceToPValue(0.0793884));
     g->SetPoint( 2 ,  900, RooStats::SignificanceToPValue(0.0129943));
     g->SetPoint( 3 ,  950, RooStats::SignificanceToPValue(0.3137820));
     g->SetPoint( 4 , 1000, RooStats::SignificanceToPValue(        0));
     g->SetPoint( 5 , 1050, RooStats::SignificanceToPValue(        0));
     g->SetPoint( 6 , 1100, RooStats::SignificanceToPValue(        0));
     g->SetPoint( 7 , 1150, RooStats::SignificanceToPValue(        0));
     g->SetPoint( 8 , 1200, RooStats::SignificanceToPValue(        0));
     g->SetPoint( 9 , 1250, RooStats::SignificanceToPValue(0.3193312));
     g->SetPoint(10 , 1300, RooStats::SignificanceToPValue(0.6791961));
     g->SetPoint(11 , 1350, RooStats::SignificanceToPValue(0.3342772));
     g->SetPoint(12 , 1400, RooStats::SignificanceToPValue(        0));
     g->SetPoint(13 , 1450, RooStats::SignificanceToPValue(        0));
     g->SetPoint(14 , 1500, RooStats::SignificanceToPValue(        0));
     g->SetPoint(15 , 1550, RooStats::SignificanceToPValue(        0));
     g->SetPoint(16 , 1600, RooStats::SignificanceToPValue(        0));
     g->SetPoint(17 , 1650, RooStats::SignificanceToPValue(        0));
     g->SetPoint(18 , 1700, RooStats::SignificanceToPValue(        0));
     g->SetPoint(19 , 1750, RooStats::SignificanceToPValue(        0));
     g->SetPoint(20 , 1800, RooStats::SignificanceToPValue(        0));
     g->SetPoint(21 , 1850, RooStats::SignificanceToPValue(        0));
     g->SetPoint(22 , 1900, RooStats::SignificanceToPValue(0.3540433));
     g->SetPoint(23 , 1950, RooStats::SignificanceToPValue(1.0621165));
     g->SetPoint(24 , 2000, RooStats::SignificanceToPValue(1.5288124));
     g->SetPoint(25 , 2050, RooStats::SignificanceToPValue(1.7860747));
     g->SetPoint(26 , 2100, RooStats::SignificanceToPValue(1.9195753));
     g->SetPoint(27 , 2150, RooStats::SignificanceToPValue(2.0153963));
     g->SetPoint(28 , 2200, RooStats::SignificanceToPValue(2.0551673));
     g->SetPoint(29 , 2250, RooStats::SignificanceToPValue(2.0020122));
     g->SetPoint(30 , 2300, RooStats::SignificanceToPValue(1.8834584));
     g->SetPoint(31 , 2350, RooStats::SignificanceToPValue(1.6989682));
     g->SetPoint(32 , 2400, RooStats::SignificanceToPValue(1.4154918));
     g->SetPoint(33 , 2450, RooStats::SignificanceToPValue(1.0144726));
     g->SetPoint(34 , 2500, RooStats::SignificanceToPValue(0.4691268));
     g->SetPoint(35 , 2550, RooStats::SignificanceToPValue(        0));
     g->SetPoint(36 , 2600, RooStats::SignificanceToPValue(        0));
     g->SetPoint(37 , 2650, RooStats::SignificanceToPValue(        0));
     g->SetPoint(38 , 2700, RooStats::SignificanceToPValue(        0));
     g->SetPoint(39 , 2750, RooStats::SignificanceToPValue(        0));
     g->SetPoint(40 , 2800, RooStats::SignificanceToPValue(        0));
     g->SetPoint(41 , 2850, RooStats::SignificanceToPValue(        0));
     g->SetPoint(42 , 2900, RooStats::SignificanceToPValue(        0));
     g->SetPoint(43 , 2950, RooStats::SignificanceToPValue(        0));
     g->SetPoint(44 , 3000, RooStats::SignificanceToPValue(        0));
     g->SetPoint(45 , 3050, RooStats::SignificanceToPValue(        0));
     g->SetPoint(46 , 3100, RooStats::SignificanceToPValue(        0));
     g->SetPoint(47 , 3150, RooStats::SignificanceToPValue(        0));
     g->SetPoint(48 , 3200, RooStats::SignificanceToPValue(        0));
     g->SetPoint(49 , 3250, RooStats::SignificanceToPValue(        0));
     g->SetPoint(50 , 3300, RooStats::SignificanceToPValue(        0));
     g->SetPoint(51 , 3350, RooStats::SignificanceToPValue(        0));
     g->SetPoint(52 , 3400, RooStats::SignificanceToPValue(        0));
     g->SetPoint(53 , 3450, RooStats::SignificanceToPValue(        0));
     g->SetPoint(54 , 3500, RooStats::SignificanceToPValue(        0));

     
     
     g->Draw("APL");
   }
     
    
    
    
    
    
    
    
    
    
    

           0 *         0 *         0 *      1000 *         1 *         0 *    123456 *         0 * 1.5543333 * 1.5569564 *        -1 *
           1 *         0 *         0 *      1050 *         1 *         0 *    123456 *         0 * 1.6186666 * 1.6214891 *        -1 *
           2 *         0 *         0 *      1100 *         1 *         0 *    123456 *         0 * 1.6643333 * 1.6673792 *        -1 *
           3 *         0 *         0 *      1150 *         1 *         0 *    123456 *         0 * 1.6879999 * 1.6923084 *        -1 *
           4 *         0 *         0 *      1200 *         1 *         0 *    123456 *         0 * 1.6495000 * 1.6513305 *        -1 *
           5 * 0.3193312 *         0 *      1250 *         1 *         0 *    123456 *         0 * 1.6805000 * 1.6821141 *        -1 *
           6 * 0.6791961 *         0 *      1300 *         1 *         0 *    123456 *         0 * 1.6260000 * 1.6273967 *        -1 *
           7 * 0.3342772 *         0 *      1350 *         1 *         0 *    123456 *         0 * 1.6898332 * 1.6913104 *        -1 *
           8 *         0 *         0 *      1400 *         1 *         0 *    123456 *         0 * 1.6208332 * 1.6224399 *        -1 *
           9 *         0 *         0 *      1450 *         1 *         0 *    123456 *         0 * 1.5913333 * 1.5929012 *        -1 *
          10 *         0 *         0 *      1500 *         1 *         0 *    123456 *         0 * 1.5163333 * 1.5177618 *        -1 *
          11 *         0 *         0 *      1550 *         1 *         0 *    123456 *         0 * 1.5505000 * 1.5519725 *        -1 *
          12 *         0 *         0 *      1600 *         1 *         0 *    123456 *         0 * 1.5043333 *  1.505759 *        -1 *
          13 *         0 *         0 *      1650 *         1 *         0 *    123456 *         0 * 1.4591666 * 1.4605387 *        -1 *
          14 *         0 *         0 *      1700 *         1 *         0 *    123456 *         0 * 1.5086666 * 1.5099616 *        -1 *
          15 *         0 *         0 *      1750 *         1 *         0 *    123456 *         0 * 1.6909999 * 1.7575948 *        -1 *
          16 *         0 *         0 *      1800 *         1 *         0 *    123456 *         0 * 1.5909999 * 1.5928876 *        -1 *
          17 *         0 *         0 *      1850 *         1 *         0 *    123456 *         0 * 1.5838333 * 1.5877807 *        -1 *
          18 * 0.3540433 *         0 *      1900 *         1 *         0 *    123456 *         0 * 1.5471667 * 1.5484924 *        -1 *
          19 * 1.0621165 *         0 *      1950 *         1 *         0 *    123456 *         0 * 1.6173332 * 1.6188288 *        -1 *
          20 * 1.5288124 *         0 *      2000 *         1 *         0 *    123456 *         0 * 1.6403332 * 1.6416810 *        -1 *
          21 * 1.7860747 *         0 *      2050 *         1 *         0 *    123456 *         0 * 1.6605000 * 1.6616907 *        -1 *
          22 * 1.9195753 *         0 *      2100 *         1 *         0 *    123456 *         0 * 1.6561666 * 1.6575751 *        -1 *
          23 * 2.0153963 *         0 *      2150 *         1 *         0 *    123456 *         0 * 1.7130000 * 1.7143733 *        -1 *
          24 * 2.0551673 *         0 *      2200 *         1 *         0 *    123456 *         0 * 1.7473332 * 1.7509790 *        -1 *
          25 * 2.0020122 *         0 *      2250 *         1 *         0 *    123456 *         0 * 1.9081666 * 1.9098042 *        -1 *
          26 * 1.8834584 *         0 *      2300 *         1 *         0 *    123456 *         0 * 1.7928333 * 1.7945222 *        -1 *
          27 * 1.6989682 *         0 *      2350 *         1 *         0 *    123456 *         0 * 1.7923333 * 1.7940772 *        -1 *
          28 * 1.4154918 *         0 *      2400 *         1 *         0 *    123456 *         0 * 1.8488333 * 1.8503024 *        -1 *
          29 * 1.0144726 *         0 *      2450 *         1 *         0 *    123456 *         0 * 1.8183333 * 1.8198983 *        -1 *
          30 * 0.4691268 *         0 *      2500 *         1 *         0 *    123456 *         0 * 1.7525000 * 1.7541691 *        -1 *
          31 *         0 *         0 *      2550 *         1 *         0 *    123456 *         0 * 1.6468333 * 1.6483798 *        -1 *
          32 *         0 *         0 *      2600 *         1 *         0 *    123456 *         0 * 1.5968333 * 1.5983787 *        -1 *
          33 *         0 *         0 *      2650 *         1 *         0 *    123456 *         0 * 1.6325000 * 1.6341334 *        -1 *
          34 *         0 *         0 *      2700 *         1 *         0 *    123456 *         0 * 1.6074999 * 1.6104178 *        -1 *
          35 *         0 *         0 *      2750 *         1 *         0 *    123456 *         0 * 1.6163333 * 1.6191442 *        -1 *
          36 *         0 *         0 *      2800 *         1 *         0 *    123456 *         0 * 1.6095000 * 1.6109803 *        -1 *
          37 *         0 *         0 *      2850 *         1 *         0 *    123456 *         0 * 1.6073333 * 1.6088748 *        -1 *
          38 *         0 *         0 *      2900 *         1 *         0 *    123456 *         0 * 1.6015000 * 1.6029496 *        -1 *
          39 *         0 *         0 *      2950 *         1 *         0 *    123456 *         0 * 1.4939999 * 1.4952561 *        -1 *
          40 *         0 *         0 *      3000 *         1 *         0 *    123456 *         0 * 1.5998333 * 1.6013703 *        -1 *
          41 *         0 *         0 *      3050 *         1 *         0 *    123456 *         0 * 1.6035000 * 1.6051152 *        -1 *
          42 *         0 *         0 *      3100 *         1 *         0 *    123456 *         0 * 1.6004999 * 1.6019206 *        -1 *
          43 *         0 *         0 *      3150 *         1 *         0 *    123456 *         0 * 2.0799999 * 2.0819680 *        -1 *
          44 *         0 *         0 *      3200 *         1 *         0 *    123456 *         0 * 1.5954999 * 1.5969563 *        -1 *
          45 *         0 *         0 *      3250 *         1 *         0 *    123456 *         0 * 1.6286666 * 1.6324663 *        -1 *
          46 *         0 *         0 *      3300 *         1 *         0 *    123456 *         0 * 1.6061667 * 1.6086767 *        -1 *
          47 *         0 *         0 *      3350 *         1 *         0 *    123456 *         0 * 1.6054999 * 1.6071476 *        -1 *
          48 *         0 *         0 *      3400 *         1 *         0 *    123456 *         0 * 1.5980000 * 1.5997926 *        -1 *
          49 *         0 *         0 *      3450 *         1 *         0 *    123456 *         0 * 1.6388332 * 1.6404980 *        -1 *
          50 *         0 *         0 *      3500 *         1 *         0 *    123456 *         0 * 1.5933333 * 1.5948863 *        -1 *
          51 *         0 *         0 *       800 *         1 *         0 *    123456 *         0 * 1.5071666 * 1.5131355 *        -1 *
          52 * 0.0793884 *         0 *       850 *         1 *         0 *    123456 *         0 * 1.5008333 * 1.5089675 *        -1 *
          53 * 0.0129943 *         0 *       900 *         1 *         0 *    123456 *         0 * 1.9084999 * 1.9116302 *        -1 *
          54 * 0.3137820 *         0 *       950 *         1 *         0 *    123456 *         0 * 1.6109999 * 1.6163628 *        -1 *
    **********************************************************************************************************************************
    
  
}