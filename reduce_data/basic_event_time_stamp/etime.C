#include <iomanip>
using namespace std;

void etime(std::string fname = "prdf_files/VSCANDATA_P00-0000359711-0000.PRDFF" ){
  gSystem->Load("libEvent.so");

  // Eventiterator *it = new fileEventiterator("/common/a1/eventdata/EVENTDATA_P00-0000410840-0078.PRDFF");
  // Eventiterator *it = new fileEventiterator("/common/p0/TestPRDFF/Run11pp500/EVENTDATA_P00-0000338145-0003.PRDFF");
  Eventiterator *it = new fileEventiterator(fname.c_str());

  Event *e;

  int i;

  int t_baseline;

  int atptime[64];
  int atpcount[64];

  for ( i = 0; i < 64; i++) 
  {
    atptime[i] = 0;
    atpcount[i] = 0;
  }


  int count = 0;

  e = it->getNextEvent();
  while (e  &&  count++ < 20000) 
  {
    //      cout << e->getFrameEntry("FRAMESOURCEID") << "  " << e->getTime() << endl;
    int atpnr = e->getFrameEntry("FRAMESOURCEID");
    if ( atpnr < 64 && atpcount[atpnr] == 0 )  
    {
      atptime[ atpnr ] +=  e->getTime() ;
      atpcount[atpnr]++;
    }
    // if( !(atpnr < 64) )
    // {
    //   cout << "ATP NUMBER " << atpnr << endl;
    // }
    delete e;
    e = it->getNextEvent();
  }

  // i is the ATP number
  for ( i = 0; i < 64; i++) 
  {

    if  (atpcount[i] ) 
    {
      cout << setw(3) <<  i << setw(12) << atptime[i] << endl;
    }
  }
  delete it;
}
