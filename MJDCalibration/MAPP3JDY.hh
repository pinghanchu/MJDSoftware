const Int_t kChannels = 58;

Int_t GoodBad0[kChannels] = {
  1,1,1,1,1,1,1,1,
  0,0,1,1,1,1,0,0,
  0,0,0,0,0,0,1,1,
  1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,0,0,1,1,
  0,0,1,1,1,1,0,0,
  1,1,1,1,1,1,0,0
};

Int_t OrtecBeGe0[kChannels] = {
  1,1,1,1,1,1,1,1,
  0,0,1,1,1,1,1,1,
  0,0,1,1,1,1,1,1,
  0,0,0,0,0,0,0,0,0,0,
  0,0,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  0,0,1,1,1,1,1,1
};

Double_t Pulser0[kChannels] = {
  0, 0, 4757.5, 1410.5, 4620.5, 1371.5, 4564.5, 1345.5,
  0, 0, 2100.5, 627.5,  431.5,  127.5,  0, 0,
  0, 0, 0, 0, 2095.5, 617.5, 2117.5, 630.5,
  4674.5, 1372.5, 405.5, 120.5, 0, 0, 2139.5, 639.5, 4574.5, 1376.5,
  4473.5, 1343.5, 2060.5, 306.5, 2099.5, 615.5, 2106.5, 626.5,
  4564.5, 1380.5, 4547.5, 1348.5, 1078.5, 316.5, 0, 0,
  4472.5, 1326.5, 4491.5, 1324.5, 4555.5, 1359.5, 0, 0
};


Double_t PulserCal0[kChannels] = {
  0., 0., 1799.5, 1798.5, 1843.5, 1843.5, 1856.5, 1855.5,
  0, 0, 844.5, 843.5, 171.5, 172.5, 0, 0,
  0, 0, 0, 0, 845.5, 844.5, 862.5, 861.5,
  1784.5, 1784.5, 154.5, 154.5, 0., 0., 813.5, 812.5, 1768.5, 1768.5,
  1829.5, 1830.5, 824.5, 822.5, 850.5, 849.5, 847.5, 847.5,
  1869.5, 1868.5, 1854.5, 1854.5, 0., 0., 0, 0,
  1847.5, 1847.5, 1837.5, 1836.5, 1834.5, 1834.5, 0.,0.
};

Double_t Pulser1[kChannels] = {
  109.5, 33.5, 1691.5, 501.5, 1600.5, 474.5, 1571.5, 463.5,
  0, 0, 1592.5, 475.5, 1599.5, 471.5, 0, 0,
  0, 0, 0, 0, 1609.5, 475.5, 1568.5, 467.5,
  1675.5, 492.5, 439.5, 131.5, 0, 0, 1671.5, 499.5, 1647.5, 495.5,
  1564.5, 469.5, 1595.5, 237.5, 1581.5, 465.5, 1587.5, 472.5,
  1564.5, 473.5, 1569.5, 465.5, 532.5, 155.5, 0, 0,
  1548.5, 459.5, 1563.5, 461.5, 1589.5, 474.5, 0, 0
};

Double_t PulserCal1[kChannels] = {
  44.5, 44.5, 639.5, 639.5, 638.5, 638.5, 639.5, 639.5,
  0, 0, 639.5, 640.5, 638.5, 637.5, 0., 0.,
  0, 0, 0, 0, 648.5, 649.5, 638.5, 638.5,
  639.5, 639.5, 167.5, 167.5, 0, 0, 635.5, 635.5, 636.5, 636.5,
  639.5, 640.5, 638.5, 637.5, 640.5, 641.5, 639.5, 639.5,
  640.5, 641.5, 639.5, 640.5, 679.5, 681.5, 0, 0,
  639.5, 640.5, 640.5, 641.5, 639.5, 639.5, 0, 0,
};
