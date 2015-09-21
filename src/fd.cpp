#include <iostream>
#include <cassert>
#include "coord.hpp"
#include "fd.hpp"
#include <mpi.h>

using namespace std;

fd_type::fd_type(const fd_type& copyfd_type) {
    
    // copy constructor
    
	h0 = copyfd_type.get_h0();
    sbporder = copyfd_type.get_sbporder();
    
    fdcoeff = new double* [2*sbporder-1];
	disscoeff = new double* [2*sbporder-1];
	
	for (int i=0; i<2*sbporder-1; i++) {
		fdcoeff[i] = new double [3*(sbporder-1)];
		disscoeff[i] = new double [3*(sbporder-1)];
	}
    
    for (int i=0; i<2*sbporder-1; i++) {
        for (int j=0; j<3*(sbporder-1); j++) {
            fdcoeff[i][j] = copyfd_type.fdcoeff[i][j];
            disscoeff[i][j] = copyfd_type.disscoeff[i][j];
        }
    }
    
}

fd_type::fd_type(int order) {
    
    // constructor to initialize coefficients based on sbporder
    
    assert(order == 2 || order == 3 || order == 4);
	
	sbporder = order;
	
	// allocate memory for FD coefficients
	
	fdcoeff = new double* [2*sbporder-1];
	disscoeff = new double* [2*sbporder-1];
	
	for (int i=0; i<2*sbporder-1; i++) {
		fdcoeff[i] = new double [3*(sbporder-1)];
		disscoeff[i] = new double [3*(sbporder-1)];
	}
	
	switch (sbporder) {
		case 2:
			// interior is 3 pt stencil
			fdcoeff[0][0] = -0.5; fdcoeff[0][1] = 0.; fdcoeff[0][2] = 0.5;
			disscoeff[0][0] = 0.5; disscoeff[0][1] = -1.; disscoeff[0][2] = 0.5;
			// boundary is 2 pts, 3 pt stencil
			fdcoeff[1][0] = -1.; fdcoeff[1][1] = 1.; fdcoeff[1][2] = 0.;
			fdcoeff[2][0] = -0.5; fdcoeff[2][1] = 0.; fdcoeff[2][2] = 0.5;
			disscoeff[1][0] = -1.; disscoeff[1][1] = 1.; disscoeff[1][2] = 0.;
			disscoeff[2][0] = 0.5; disscoeff[2][1] = -1.; disscoeff[2][2] = 0.5;
			// boundary norm
			h0 = 0.5;
			break;
		case 3:
			// interior is 5 pt stencil
			fdcoeff[0][0] = 1./12.; fdcoeff[0][1] = -2./3.; fdcoeff[0][2] = 0.; fdcoeff[0][3] = 2./3.; fdcoeff[0][4] = -1./12.;
			disscoeff[0][0] = -1./12.; disscoeff[0][1] = 1./3.; disscoeff[0][2] = -1./2.; disscoeff[0][3] = 1./3.; disscoeff[0][4] = -1./12.;
			// boundary is 4 pts, 6 pt stencil
			fdcoeff[1][0] = -24./17.; fdcoeff[1][1] = 59./34.; fdcoeff[1][2] = -4./17.; fdcoeff[1][3] = -3./34.; fdcoeff[1][4] = 0.; fdcoeff[1][5] = 0.;
			fdcoeff[2][0] = -0.5; fdcoeff[2][1] = 0.; fdcoeff[2][2] = 0.5; fdcoeff[2][3] = 0.; fdcoeff[2][4] = 0.; fdcoeff[2][5] = 0.;
			fdcoeff[3][0] = 4./43.; fdcoeff[3][1] = -59./86.; fdcoeff[3][2] = 0.; fdcoeff[3][3] = 59./86.; fdcoeff[3][4] = -4./43.; fdcoeff[3][5] = 0.;
			fdcoeff[4][0] = 3./98.; fdcoeff[4][1] = 0.; fdcoeff[4][2] = -59./98.; fdcoeff[4][3] = 0.; fdcoeff[4][4] = 32./49.; fdcoeff[4][5] = -4./49.;
			disscoeff[1][0] = -8./17.; disscoeff[1][1] = 16./17.; disscoeff[1][2] = -8./17.; disscoeff[1][3] = 0.; disscoeff[1][4] = 0.; disscoeff[1][5] = 0.;
			disscoeff[2][0] = 16./59.; disscoeff[2][1] = -36./59.; disscoeff[2][2] = 24./59.; disscoeff[2][3] = -4./59.; disscoeff[2][4] = 0.; disscoeff[2][5] = 0.;
			disscoeff[3][0] = -8./43.; disscoeff[3][1] = 24./43.; disscoeff[3][2] = -28./43.; disscoeff[3][3] = 16./43.; disscoeff[3][4] = -4./43; disscoeff[3][5] = 0.;
			disscoeff[4][0] = 0.; disscoeff[4][1] = -4./49.; disscoeff[4][2] = 16./49.; disscoeff[4][3] = -24./49.; disscoeff[4][4] = 16./49.; disscoeff[4][5] = -4./49.;
			// boundary norm
			h0 = 17./48.;
			break;
		case 4:
			// interior is 7 pt stencil
			fdcoeff[0][0] = -1./60.; fdcoeff[0][1] = 3./20.; fdcoeff[0][2] = -3./4.; fdcoeff[0][3] = 0.; fdcoeff[0][4] = 3./4.; fdcoeff[0][5] = -3./20.; fdcoeff[0][6] = 1./60.;
			disscoeff[0][0] = 1./60.; disscoeff[0][1] = -1./10.; disscoeff[0][2] = 1./4.; disscoeff[0][3] = -1./3.; disscoeff[0][4] = 1./4.; disscoeff[0][5] = -1./10.; disscoeff[0][6] = 1./60.;
			// boundary is 6 pts, 9 pt stencil
            
            // version from Gustaffson, 2008

/*            fdcoeff[1][0] = -21600./13649.; fdcoeff[1][1] = 104009./54596.; fdcoeff[1][2] = 30443./81894.; fdcoeff[1][3] = -33311./27298.; fdcoeff[1][4] = 16863./27298.; fdcoeff[1][5] = -15025./163788.;   fdcoeff[1][6] = 0.; fdcoeff[1][7] = 0.; fdcoeff[1][8] = 0.;
            fdcoeff[2][0] = -104009./240260.; fdcoeff[2][1] = 0.; fdcoeff[2][2] = -311./72078.; fdcoeff[2][3] = 20229./24026.; fdcoeff[2][4] = -24337./48052.; fdcoeff[2][5] = 36661./360390.; fdcoeff[2][6] = 0.; fdcoeff[2][7] = 0.; fdcoeff[2][8] = 0.;
            fdcoeff[3][0] = -30443./162660.; fdcoeff[3][1] = 311./32532.; fdcoeff[3][2] = 0.; fdcoeff[3][3] = -11155./16266.; fdcoeff[3][4] = 41287./32532.; fdcoeff[3][5] = -21999./54220.; fdcoeff[3][6] = 0.; fdcoeff[3][7] = 0.; fdcoeff[3][8] = 0.;
            fdcoeff[4][0] = 33311./107180.; fdcoeff[4][1] = -20229./21436.; fdcoeff[4][2] = 485./1398.; fdcoeff[4][3] = 0.; fdcoeff[4][4] = 4147./21436.; fdcoeff[4][5] = 25427./321540.; fdcoeff[4][6] = 72./5359.; fdcoeff[4][7] = 0.; fdcoeff[4][8] = 0.;
            fdcoeff[5][0] = -16863./78770.; fdcoeff[5][1] = 24337./31508.; fdcoeff[5][2] = -41287./47262.; fdcoeff[5][3] = -4147./15754.; fdcoeff[5][4] = 0.; fdcoeff[5][5] = 342523./472620.; fdcoeff[5][6] = -1296./7877.; fdcoeff[5][7] = 144./7877.; fdcoeff[5][8] = 0.;
            fdcoeff[6][0] = 15025./525612.; fdcoeff[6][1] = -36661./262806.; fdcoeff[6][2] = 21999./87602.; fdcoeff[6][3] = -25427./262806.; fdcoeff[6][4] = -342523./525612.; fdcoeff[6][5] = 0.; fdcoeff[6][6] = 32400./43801.; fdcoeff[6][7] = -6480./43801.; fdcoeff[6][8] = 720./43801.;*/
            
            // version from fdmap code (unsure of source) with dissipation operator
            
/*			fdcoeff[1][0] = -1.582533518939116; fdcoeff[1][1] = 2.033378678700676; fdcoeff[1][2] = -0.141512858744873; fdcoeff[1][3] = -0.450398306578272; fdcoeff[1][4] = 0.104488069284042; fdcoeff[1][5] = 0.036577936277544; fdcoeff[1][6] = 0.; fdcoeff[1][7] = 0.; fdcoeff[1][8] = 0.;
			fdcoeff[2][0] = -0.462059195631158; fdcoeff[2][1] = 0.; fdcoeff[2][2] = 0.287258622978251; fdcoeff[2][3] = 0.258816087376832; fdcoeff[2][4] = -0.069112065532624; fdcoeff[2][5] = -0.014903449191300; fdcoeff[2][6] = 0.; fdcoeff[2][7] = 0.; fdcoeff[2][8] = 0.;
			fdcoeff[3][0] = 0.071247104721830; fdcoeff[3][1] = -0.636451095137907; fdcoeff[3][2] = 0.; fdcoeff[3][3] = 0.606235523609147; fdcoeff[3][4] = -0.022902190275815; fdcoeff[3][5] = -0.018129342917256; fdcoeff[3][6] = 0.; fdcoeff[3][7] = 0.; fdcoeff[3][8] = 0.;
			fdcoeff[4][0] = 0.114713313798970; fdcoeff[4][1] = -0.290087484386815; fdcoeff[4][2] = -0.306681191361148; fdcoeff[4][3] = 0.; fdcoeff[4][4] = 0.520262285050482; fdcoeff[4][5] = -0.051642265516119; fdcoeff[4][6] = 0.013435342414630; fdcoeff[4][7] = 0.; fdcoeff[4][8] = 0.;
			fdcoeff[5][0] = -0.036210680656541; fdcoeff[5][1] = 0.105400944933782; fdcoeff[5][2] = 0.015764336127392; fdcoeff[5][3] = -0.707905442575989; fdcoeff[5][4] = 0.; fdcoeff[5][5] = 0.769199413962647; fdcoeff[5][6] = -0.164529643265203; fdcoeff[5][7] = 0.018281071473911; fdcoeff[5][8] = 0.;
			fdcoeff[6][0] = -0.011398193015050; fdcoeff[6][1] = 0.020437334208704; fdcoeff[6][2] = 0.011220896474665; fdcoeff[6][3] = 0.063183694641876; fdcoeff[6][4] = -0.691649024426814; fdcoeff[6][5] = 0.; fdcoeff[6][6] = 0.739709139060752; fdcoeff[6][7] = -0.147941827812150; fdcoeff[6][8] = 0.016437980868017;
			disscoeff[1][0] = -0.105502234595941; disscoeff[1][1] = 0.316506703787823; disscoeff[1][2] = -0.316506703787823; disscoeff[1][3] = 0.105502234595941; disscoeff[1][4] = 0.; disscoeff[1][5] = 0.; disscoeff[1][6] = 0.; disscoeff[1][7] = 0.; disscoeff[1][8] = 0.;
			disscoeff[2][0] = 0.071922084408557; disscoeff[2][1] = -0.227753267293765; disscoeff[2][2] = 0.251727295429951; disscoeff[2][3] = -0.107883126612836; disscoeff[2][4] = 0.011987014068093; disscoeff[2][5] = 0.; disscoeff[2][6] = 0.; disscoeff[2][7] = 0.; disscoeff[2][8] = 0.;
			disscoeff[3][0] = -0.159350793065290; disscoeff[3][1] = 0.557727775728513; disscoeff[3][2] = -0.743637034304685; disscoeff[3][3] = 0.478052379195869; disscoeff[3][4] = -0.159350793065290; disscoeff[3][5] = 0.026558465510882; disscoeff[3][6] = 0.; disscoeff[3][7] = 0.; disscoeff[3][8] = 0.;
			disscoeff[4][0] = 0.026870684829259; disscoeff[4][1] = -0.120918081731666; disscoeff[4][2] = 0.241836163463333; disscoeff[4][3] = -0.282142190707221; disscoeff[4][4] = 0.201530136219444; disscoeff[4][5] = -0.080612054487778; disscoeff[4][6] = 0.013435342414630; disscoeff[4][7] = 0.; disscoeff[4][8] = 0.;
			disscoeff[5][0] = 0.; disscoeff[5][1] = 0.018281071473911; disscoeff[5][2] = -0.109686428843468; disscoeff[5][3] = 0.274216072108671; disscoeff[5][4] = -0.365621429478228; disscoeff[5][5] = 0.274216072108671; disscoeff[5][6] = -0.109686428843468; disscoeff[5][7] = 0.018281071473911; disscoeff[5][8] = 0.;
            disscoeff[6][0] = 0.; disscoeff[6][1] = 0.; disscoeff[6][2] = 0.016437980868017; disscoeff[6][3] = -0.098627885208100; disscoeff[6][4] = 0.246569713020251; disscoeff[6][5] = -0.328759617360334; disscoeff[6][6] = 0.246569713020251; disscoeff[6][7] = -0.098627885208100; disscoeff[6][8] = 0.016437980868017;*/
            
            // minimum width version from Strand (1994)
            
			fdcoeff[1][0] = -21600./13649.; fdcoeff[1][1] = 81763./40947.; fdcoeff[1][2] = 131./27298.; fdcoeff[1][3] = -9143./13649.; fdcoeff[1][4] = 20539./81894.; fdcoeff[1][5] = 0.; fdcoeff[1][6] = 0.; fdcoeff[1][7] = 0.; fdcoeff[1][8] = 0.;
            fdcoeff[2][0] = -81763./180195.; fdcoeff[2][1] = 0.; fdcoeff[2][2] = 7357./36039.; fdcoeff[2][3] = 30637./72078.; fdcoeff[2][4] = -2328./12013.; fdcoeff[2][5] = 6611./360390.; fdcoeff[2][6] = 0.; fdcoeff[2][7] = 0.; fdcoeff[2][8] = 0.;
            fdcoeff[3][0] = -131./54220.; fdcoeff[3][1] = -7357./16266.; fdcoeff[3][2] = 0.; fdcoeff[3][3] = 645./2711.; fdcoeff[3][4] = 11237./32532.; fdcoeff[3][5] = -3487./27110.; fdcoeff[3][6] = 0.; fdcoeff[3][7] = 0.; fdcoeff[3][8] = 0.;
            fdcoeff[4][0] = 9143./53590.; fdcoeff[4][1] = -30637./64308.; fdcoeff[4][2] = -645./5359.; fdcoeff[4][3] = 0.; fdcoeff[4][4] = 13733./32154.; fdcoeff[4][5] = -67./4660.; fdcoeff[4][6] = 72./5359.; fdcoeff[4][7] = 0.; fdcoeff[4][8] = 0.;
            fdcoeff[5][0] = -20539./236310.; fdcoeff[5][1] = 2328./7877.; fdcoeff[5][2] = -11237./47262.; fdcoeff[5][3] = -13733./23631.; fdcoeff[5][4] = 0.; fdcoeff[5][5] = 89387./118155.; fdcoeff[5][6] = -1296./7877.; fdcoeff[5][7] = 144./7877.; fdcoeff[5][8] = 0.;
            fdcoeff[6][0] = 0.; fdcoeff[6][1] = -6611./262806.; fdcoeff[6][2] = 3487./43801.; fdcoeff[6][3] = 1541./87602.; fdcoeff[6][4] = -89387./131403.; fdcoeff[6][5] = 0.; fdcoeff[6][6] = 32400./43801.; fdcoeff[6][7] = -6480./43801.; fdcoeff[6][8] = 720./43801.;
            
			// boundary norm
			h0 = 13649./43200.;
			break;
        case 5:
            fdcoeff[0][0] = 1./280.; fdcoeff[0][1] = -4./105.; fdcoeff[0][2] = 1./5.; fdcoeff[0][3] = -4./5.; fdcoeff[0][4] = 0.; fdcoeff[0][5] = 4./5.; fdcoeff[0][6] = -1./5.; fdcoeff[0][7] = 4./105.; fdcoeff[0][8] = -1./280.;
            fdcoeff[1][0] = -2540160./1498139.; fdcoeff[1][1] = 699846290793./311403172540.; fdcoeff[1][2] = -10531586157./311403172540.; fdcoeff[1][3] = -145951651079./186841903524.; fdcoeff[1][4] = 398124597./15570158627.; fdcoeff[1][5] = 39152001./113858564.; fdcoeff[1][6] = -80631892229./934209517620.; fdcoeff[1][7] =  -6230212503./311403172540.; fdcoeff[1][8] =  0.; fdcoeff[1][9] =  0.; fdcoeff[1][10] =  0.; fdcoeff[1][11] =  0.;
            fdcoeff[2][0] =  -24132630717./55557028660.; fdcoeff[2][1] =  0.; fdcoeff[2][2] =  2113176981./23016483302.; fdcoeff[2][3] =  5686186719./11508241651.; fdcoeff[2][4] =  -3408473341./138098899812.; fdcoeff[2][5] =  -39291999./210388330.; fdcoeff[2][6] =  607046586./11508241651.; fdcoeff[2][7] =  3460467023./483346149342.; fdcoeff[2][8] =  0.; fdcoeff[2][9] =  0.; fdcoeff[2][10] =  0.; fdcoeff[2][11] =  0.;
            fdcoeff[3][0] =  3510528719./90623010660.; fdcoeff[3][1] =  -704392327./1294614438.; fdcoeff[3][2] =  0.; fdcoeff[3][3] =  503511235./1294614438.; fdcoeff[3][4] =  354781619./2589228876.; fdcoeff[3][5] =  407439./3944590.; fdcoeff[3][6] =  -2986842./16597621.; fdcoeff[3][7] =  169381493./3020767022.; fdcoeff[3][8] =  0.; fdcoeff[3][9] =  0.; fdcoeff[3][10] =  0.; fdcoeff[3][11] =  0.;
            fdcoeff[4][0] =  145951651079./1139279786988.; fdcoeff[4][1] =  -5686186719./13562854607.; fdcoeff[4][2] =  -1510533705./27125709214.; fdcoeff[4][3] =  0.; fdcoeff[4][4] =  6763379967./54251418428.; fdcoeff[4][5] =  13948923./49589962.; fdcoeff[4][6] =  -1603900430./40688563821.; fdcoeff[4][7] =  -3742312557./189879964498.; fdcoeff[4][8] =  0.; fdcoeff[4][9] =  0.; fdcoeff[4][10] =  0.; fdcoeff[4][11] =  0.;
            fdcoeff[5][0] =  -398124597./21790888777.; fdcoeff[5][1] =  3408473341./37355809332.; fdcoeff[5][2] =  -1064344857./12451936444.; fdcoeff[5][3] =  -6763379967./12451936444.; fdcoeff[5][4] =  0.; fdcoeff[5][5] =  763665./1198108.; fdcoeff[5][6] =  -1282435899./12451936444.; fdcoeff[5][7] =  7822226819./261490665324.; fdcoeff[5][8] =  -2592./299527.; fdcoeff[5][9] =  0.; fdcoeff[5][10] =  0.; fdcoeff[5][11] =  0.;
            fdcoeff[6][0] =  -1864381./23506116.; fdcoeff[6][1] =  13097333./58765290.; fdcoeff[6][2] =  -407439./19588430.; fdcoeff[6][3] =  -4649641./11753058.; fdcoeff[6][4] =  -254555./1237164.; fdcoeff[6][5] =  0.; fdcoeff[6][6] =  5346432./9794215.; fdcoeff[6][7] =  -923328./9794215.; fdcoeff[6][8] =  3072./103097.; fdcoeff[6][9] =  -288./103097.; fdcoeff[6][10] =  0.; fdcoeff[6][11] =  0.;
            fdcoeff[7][0] =  11518841747./417855345780.; fdcoeff[7][1] =  -607046586./6964255763.; fdcoeff[7][2] =  1192698./23768791.; fdcoeff[7][3] =  1603900430./20892767289.; fdcoeff[7][4] =  1282435899./27857023052.; fdcoeff[7][5] =  -48117888./63658645.; fdcoeff[7][6] =  0.; fdcoeff[7][7] =  301190400./366539777.; fdcoeff[7][8] =  -145152./670091.; fdcoeff[7][9] =  27648./670091.; fdcoeff[7][10] =  -2592./670091.; fdcoeff[7][11] =  0.;
            fdcoeff[8][0] =  6230212503./1065851828540.; fdcoeff[8][1] =  -3460467023./319755548562.; fdcoeff[8][2] =  -1524433437./106585182854.; fdcoeff[8][3] =  3742312557./106585182854.; fdcoeff[8][4] =  -7822226819./639511097124.; fdcoeff[8][5] =  58169664./487135205.; fdcoeff[8][6] =  -2108332800./2804873233.; fdcoeff[8][7] =  0.; fdcoeff[8][8] =  4064256./5127739.; fdcoeff[8][9] =  -1016064./5127739.; fdcoeff[8][10] =  193536./5127739.; fdcoeff[8][11] =  -18144./5127739.;
            h0 = 1498139./5080320.;
            break;
        case 6:
            fdcoeff[0][0] = -1./1260.; fdcoeff[0][1] = 5./504.; fdcoeff[0][2] = -5./84.; fdcoeff[0][3] = 5./21.; fdcoeff[0][4] = -5./6.; fdcoeff[0][5] = 0.; fdcoeff[0][6] = 5./6.; fdcoeff[0][7] = -5./21.; fdcoeff[0][8] = 5./84.; fdcoeff[0][9] = -5./504.; fdcoeff[0][10] = 1./1260.;
            fdcoeff[1][0] =  -1.7380923775745425e+00; fdcoeff[1][1] =  2.3557601935237220e+00; fdcoeff[1][2] =  -1.5328406598563976e-01; fdcoeff[1][3] =  -5.7266565770416333e-01; fdcoeff[1][4] =  -1.8308103515008173e-01; fdcoeff[1][5] =  1.8186748267946842e-01; fdcoeff[1][6] =  2.0034232582598244e-01; fdcoeff[1][7] =  2.2678007363666621e-02; fdcoeff[1][8] =  -1.1782459320459637e-01; fdcoeff[1][9] =  -3.0591175636402144e-02; fdcoeff[1][10] =  3.4890895862586133e-02; fdcoeff[1][11] =  0.0000000000000000e+00; fdcoeff[1][12] =  0.0000000000000000e+00; fdcoeff[1][13] =  0.0000000000000000e+00; fdcoeff[1][14] =  0.0000000000000000e+00; fdcoeff[1][15] =  0.0000000000000000e+00;
            fdcoeff[2][0] =  -4.3020203737210871e-01; fdcoeff[2][1] =  0.0000000000000000e+00; fdcoeff[2][2] =  1.1837297346927406e-01; fdcoeff[2][3] =  3.3928601158526644e-01; fdcoeff[2][4] =  1.3241927733034406e-01; fdcoeff[2][5] =  -8.7495003780608913e-02; fdcoeff[2][6] =  -1.1750484124279399e-01; fdcoeff[2][7] =  -1.6401912273575153e-02; fdcoeff[2][8] =  6.2537843443041474e-02; fdcoeff[2][9] =  1.7143274696828435e-02; fdcoeff[2][10] =  -1.8155585855667674e-02; fdcoeff[2][11] =  0.0000000000000000e+00; fdcoeff[2][12] =  0.0000000000000000e+00; fdcoeff[2][13] =  0.0000000000000000e+00; fdcoeff[2][14] =  0.0000000000000000e+00; fdcoeff[2][15] =  0.0000000000000000e+00;
            fdcoeff[3][0] =  3.4348531361887280e-01; fdcoeff[3][1] =  -1.4525207124434036e+00; fdcoeff[3][2] =  0.0000000000000000e+00; fdcoeff[3][3] =  2.9011513992277767e+00; fdcoeff[3][4] =  -2.2419288742360557e+00; fdcoeff[3][5] =  -5.4662873578741478e-01; fdcoeff[3][6] =  1.2908050607446131e+00; fdcoeff[3][7] =  6.1514504292452719e-02; fdcoeff[3][8] =  -4.2442625460011202e-01; fdcoeff[3][9] =  1.5579158905288801e-02; fdcoeff[3][10] =  5.2969140277981920e-02; fdcoeff[3][11] =  0.0000000000000000e+00; fdcoeff[3][12] =  0.0000000000000000e+00; fdcoeff[3][13] =  0.0000000000000000e+00; fdcoeff[3][14] =  0.0000000000000000e+00; fdcoeff[3][15] =  0.0000000000000000e+00;
            fdcoeff[4][0] =  8.6111387816878188e-02; fdcoeff[4][1] =  -2.7937273515056432e-01; fdcoeff[4][2] =  -1.9467880944770807e-01; fdcoeff[4][3] =  0.0000000000000000e+00; fdcoeff[4][4] =  2.0170150914578375e-01; fdcoeff[4][5] =  2.4269917331475005e-01; fdcoeff[4][6] =  -7.7261988327590472e-02; fdcoeff[4][7] =  5.0649247607525059e-02; fdcoeff[4][8] =  -7.4775049946661561e-03; fdcoeff[4][9] =  -4.0978487203372188e-02; fdcoeff[4][10] =  1.8608207238964152e-02; fdcoeff[4][11] =  0.0000000000000000e+00; fdcoeff[4][12] =  0.0000000000000000e+00; fdcoeff[4][13] =  0.0000000000000000e+00; fdcoeff[4][14] =  0.0000000000000000e+00; fdcoeff[4][15] =  0.0000000000000000e+00;
            fdcoeff[5][0] =  9.1509035082611684e-02; fdcoeff[5][1] =  -3.6243526359648576e-01; fdcoeff[5][2] =  5.0007055839856984e-01; fdcoeff[5][3] =  -6.7045605191055857e-01; fdcoeff[5][4] =  0.0000000000000000e+00; fdcoeff[5][5] =  -1.7807807859119628e-02; fdcoeff[5][6] =  7.5000761407401195e-01; fdcoeff[5][7] =  -2.2979723229714316e-01; fdcoeff[5][8] =  -1.2521154324370892e-01; fdcoeff[5][9] =  6.8278284106004450e-02; fdcoeff[5][10] =  -4.1575927541817690e-03; fdcoeff[5][11] =  0.0000000000000000e+00; fdcoeff[5][12] =  0.0000000000000000e+00; fdcoeff[5][13] =  0.0000000000000000e+00; fdcoeff[5][14] =  0.0000000000000000e+00; fdcoeff[5][15] =  0.0000000000000000e+00;
            fdcoeff[6][0] =  -7.5752056274147259e-02; fdcoeff[6][1] =  1.9956355926115746e-01; fdcoeff[6][2] =  1.0160630736447970e-01; fdcoeff[6][3] =  -6.7227694623145351e-01; fdcoeff[6][4] =  1.4839839882599690e-02; fdcoeff[6][5] =  0.0000000000000000e+00; fdcoeff[6][6] =  5.4091068834671807e-01; fdcoeff[6][7] =  -1.2712520372174399e-01; fdcoeff[6][8] =  -8.9292453564020990e-02; fdcoeff[6][9] =  1.6181541970619609e-01; fdcoeff[6][10] =  -5.4289154769785249e-02; fdcoeff[6][11] =  0.0000000000000000e+00; fdcoeff[6][12] =  0.0000000000000000e+00; fdcoeff[6][13] =  0.0000000000000000e+00; fdcoeff[6][14] =  0.0000000000000000e+00; fdcoeff[6][15] =  0.0000000000000000e+00;
            fdcoeff[7][0] =  -3.3838029883391296e-02; fdcoeff[7][1] =  1.0867927550524317e-01; fdcoeff[7][2] =  -9.7293058702223670e-02; fdcoeff[7][3] =  8.6783825404790446e-02; fdcoeff[7][4] =  -2.5344131542932297e-01; fdcoeff[7][5] =  -2.1934035945002228e-01; fdcoeff[7][6] =  0.0000000000000000e+00; fdcoeff[7][7] =  2.7184438867288430e-01; fdcoeff[7][8] =  1.9102691945078512e-01; fdcoeff[7][9] =  -4.8646826827046824e-02; fdcoeff[7][10] =  -6.2407959378425991e-03; fdcoeff[7][11] =  4.6597719614658163e-04; fdcoeff[7][12] =  0.0000000000000000e+00; fdcoeff[7][13] =  0.0000000000000000e+00; fdcoeff[7][14] =  0.0000000000000000e+00; fdcoeff[7][15] =  0.0000000000000000e+00;
            fdcoeff[8][0] =  -1.5567948806367624e-02; fdcoeff[8][1] =  6.1656604470023607e-02; fdcoeff[8][2] =  -1.8844858059892756e-02; fdcoeff[8][3] =  -2.3122780265804038e-01; fdcoeff[8][4] =  3.1560994521078772e-01; fdcoeff[8][5] =  2.0951677187991255e-01; fdcoeff[8][6] =  -1.1048784865195491e+00; fdcoeff[8][7] =  0.0000000000000000e+00; fdcoeff[8][8] =  1.1823059621092409e+00; fdcoeff[8][9] =  -5.3610400867086083e-01; fdcoeff[8][10] =  1.5931375952374752e-01; fdcoeff[8][11] =  -2.3673846172827626e-02; fdcoeff[8][12] =  1.8939076938262100e-03; fdcoeff[8][13] =  0.0000000000000000e+00; fdcoeff[8][14] =  0.0000000000000000e+00; fdcoeff[8][15] =  0.0000000000000000e+00;
            fdcoeff[9][0] =  2.6737701764454301e-02; fdcoeff[9][1] =  -7.7712278574126673e-02; fdcoeff[9][2] =  4.2981266272823705e-02; fdcoeff[9][3] =  1.1284579710276557e-02; fdcoeff[9][4] =  5.6847566375570611e-02; fdcoeff[9][5] =  4.8647834370398067e-02; fdcoeff[9][6] =  -2.5665536068472994e-01; fdcoeff[9][7] =  -3.9083324869946684e-01; fdcoeff[9][8] =  0.0000000000000000e+00; fdcoeff[9][9] =  6.5716944195909766e-01; fdcoeff[9][10] =  -1.5822272208022428e-01; fdcoeff[9][11] =  4.6954983762905661e-02; fdcoeff[9][12] =  -7.8258306271509429e-03; fdcoeff[9][13] =  6.2606645017207550e-04; fdcoeff[9][14] =  0.0000000000000000e+00; fdcoeff[9][15] =  0.0000000000000000e+00;
            fdcoeff[10][0] =  9.4425181052687698e-03; fdcoeff[10][1] =  -2.8976375375532045e-02; fdcoeff[10][2] =  -2.1459742428921558e-03; fdcoeff[10][3] =  8.4117843695442701e-02; fdcoeff[10][4] =  -4.2165149106440383e-02; fdcoeff[10][5] =  -1.1991463562335723e-01; fdcoeff[10][6] =  8.8902467992349743e-02; fdcoeff[10][7] =  2.4105392677971343e-01; fdcoeff[10][8] =  -8.9388344421253152e-01; fdcoeff[10][9] =  0.0000000000000000e+00; fdcoeff[10][10] =  8.6496680152924643e-01; fdcoeff[10][11] =  -2.5547312415382800e-01; fdcoeff[10][12] =  6.3868281038457000e-02; fdcoeff[10][13] =  -1.0644713506409501e-02; fdcoeff[10][14] =  8.5157708051276015e-04; fdcoeff[10][15] =  0.0000000000000000e+00;
            fdcoeff[11][0] =  -9.9625965676187218e-03; fdcoeff[11][1] =  2.8387641187789508e-02; fdcoeff[11][2] =  -6.7495090936003027e-03; fdcoeff[11][3] =  -3.5335033597892078e-02; fdcoeff[11][4] =  2.3750992019053968e-03; fdcoeff[11][5] =  3.7216380474824604e-02; fdcoeff[11][6] =  1.0550378667904333e-02; fdcoeff[11][7] =  -6.6265458456725809e-02; fdcoeff[11][8] =  1.9908619649258188e-01; fdcoeff[11][9] =  -8.0014409359906680e-01; fdcoeff[11][10] =  0.0000000000000000e+00; fdcoeff[11][11] =  8.2714572225493910e-01; fdcoeff[11][12] =  -2.3632734921569687e-01; fdcoeff[11][13] =  5.9081837303924217e-02; fdcoeff[11][14] =  -9.8469728839873684e-03; fdcoeff[11][15] =  7.8775783071898962e-04;
            h0 = 5261271563./18289152000.;
            break;
		default:
			cerr << "Error in initializing FD method\n";
			MPI_Abort(MPI_COMM_WORLD, -1);
			break;
	}

}

fd_type::~fd_type() {
	// destructor, deallocates memory for coefficients
	for (int i=0; i<2*sbporder-1; i++){
		delete[] fdcoeff[i];
		delete[] disscoeff[i];
	}
	delete[] fdcoeff;
	delete[] disscoeff;
}

fd_type& fd_type:: operator=(const fd_type& assignfd_type) {
 
    // assignment operator
    
    if (this != &assignfd_type) {
    
        h0 = assignfd_type.get_h0();
        sbporder = assignfd_type.get_sbporder();

        for (int i=0; i<2*sbporder-1; i++) {
            for (int j=0; j<3*(sbporder-1); j++) {
                fdcoeff[i][j] = assignfd_type.fdcoeff[i][j];
                disscoeff[i][j] = assignfd_type.disscoeff[i][j];
            }
        }
    }
    
    return *this;
}

int fd_type::get_sbporder() const {
    // returns sbporder
    
    return sbporder;
}

double fd_type::get_h0() const {
    // returns h0

    return h0;
}

double fd_type::cons_s(double**** f, double**** m, double*** jac, const int i, const int j, const int k, const coord c, const int dir, const int index, const int ndim, const int mode) const {
    // returns conservative finite difference of the given field for the indices given
    // in the specified direction
    
/*    assert(index == 0 || index == 1 || index == 2);
    assert(ndim == 2 || ndim == 3);
    assert(mode == 2 || mode == 3);
    assert(dir == 0 || dir == 1 || dir == 2);
    assert(i >=0 && i < c.get_nx_tot(0));
    assert(j >=0 && j < c.get_nx_tot(1));
    assert(k >=0 && k < c.get_nx_tot(2));*/
    
    double fdiff = 0.;
    int n, n_loc, m_loc, m_ghost;
    int findex[3];
    
    switch (ndim) {
    case 3:
        switch (index) {
        case 0:
            findex[0] = 0;
            findex[1] = 1;
            findex[2] = 2;
            break;
        case 1:
            findex[0] = 1;
            findex[1] = 3;
            findex[2] = 4;
            break;
        case 2:
            findex[0] = 2;
            findex[1] = 4;
            findex[2] = 5;
        }
        break;
    case 2:
        switch (mode) {
        case 2:
            findex[0] = 0;
            findex[1] = 1;
            break;
        case 3:
            switch (index) {
            case 0:
                findex[0] = 0;
                findex[1] = 1;
                break;
            case 1:
                findex[0] = 1;
                findex[1] = 2;
            }
        }
    }
    
    switch (dir) {
    case 0:
        n = c.get_nx(0);
        n_loc = c.get_nx_loc(0);
        m_loc = c.get_xm_loc(0)-c.get_xm(0);
        m_ghost = c.get_xm_ghost(0);
        if (i<2*(sbporder-1) && m_loc == 0) {
            for (int ii=0; ii<3*(sbporder-1); ii++) {
                for (int iii=0; iii<ndim; iii++) {
                    fdiff += (fdcoeff[i+1][ii]*jac[ii][j][k]*m[iii][ii][j][k]*
                              f[findex[iii]][ii][j][k]);
                }
            }
        } else if (i>n_loc+m_ghost-1-2*(sbporder-1) && m_loc+n_loc == n) {
            for (int ii=0; ii<3*(sbporder-1); ii++) {
                for (int iii=0; iii<ndim; iii++) {
                    fdiff -= (fdcoeff[n_loc+m_ghost-i][ii]*jac[n_loc+m_ghost-1-ii][j][k]*
                              m[iii][n_loc+m_ghost-1-ii][j][k]*
                              f[findex[iii]][n_loc+m_ghost-1-ii][j][k]);
                }
            }
        } else {
            for (int ii=0; ii<2*(sbporder-1)+1; ii++) {
                for (int iii=0; iii<ndim; iii++) {
                    fdiff += (fdcoeff[0][ii]*jac[i-sbporder+1+ii][j][k]*
                              m[iii][i-sbporder+1+ii][j][k]*
                              f[findex[iii]][i-sbporder+1+ii][j][k]);
                }
            }
        }
        break;
    case 1:
        n = c.get_nx(1);
        n_loc = c.get_nx_loc(1);
        m_loc = c.get_xm_loc(1)-c.get_xm(1);
        m_ghost = c.get_xm_ghost(1);
        if (j<2*(sbporder-1) && m_loc == 0) {
            for (int jj=0; jj<3*(sbporder-1); jj++) {
                for (int jjj=0; jjj<ndim; jjj++) {
                    fdiff += (fdcoeff[j+1][jj]*jac[i][jj][k]*m[jjj][i][jj][k]*
                              f[findex[jjj]][i][jj][k]);
                }
            }
        } else if (j>n_loc+m_ghost-1-2*(sbporder-1) && m_loc+n_loc == n) {
            for (int jj=0; jj<3*(sbporder-1); jj++) {
                for (int jjj=0; jjj<ndim; jjj++) {
                    fdiff -= (fdcoeff[n_loc+m_ghost-j][jj]*jac[i][n_loc+m_ghost-1-jj][k]*
                              m[jjj][i][n_loc+m_ghost-1-jj][k]*
                              f[findex[jjj]][i][n_loc+m_ghost-1-jj][k]);
                }
            }
        } else {
            for (int jj=0; jj<2*(sbporder-1)+1; jj++) {
                for (int jjj=0; jjj<ndim; jjj++) {
                    fdiff += (fdcoeff[0][jj]*jac[i][j-sbporder+1+jj][k]*
                              m[jjj][i][j-sbporder+1+jj][k]*
                              f[findex[jjj]][i][j-sbporder+1+jj][k]);
                }
            }
        }
        break;
    case 2:
        n = c.get_nx(2);
        n_loc = c.get_nx_loc(2);
        m_loc = c.get_xm_loc(2)-c.get_xm(2);
        m_ghost = c.get_xm_ghost(2);
        if (k<2*(sbporder-1) && m_loc == 0) {
            for (int kk=0; kk<3*(sbporder-1); kk++) {
                for (int kkk=0; kkk<ndim; kkk++) {
                    fdiff += (fdcoeff[k+1][kk]*jac[i][j][kk]*m[kkk][i][j][kk]*
                              f[findex[kkk]][i][j][kk]);
                }
            }
        } else if (k>n_loc+m_ghost-1-2*(sbporder-1) && m_loc+n_loc == n) {
            for (int kk=0; kk<3*(sbporder-1); kk++) {
                for (int kkk=0; kkk<ndim; kkk++) {
                    fdiff -= (fdcoeff[n_loc+m_ghost-k][kk]*jac[i][j][n_loc+m_ghost-1-kk]*
                              m[kkk][i][j][n_loc+m_ghost-1-kk]*
                              f[findex[kkk]][i][j][n_loc+m_ghost-1-kk]);
                }
            }
        } else {
            for (int kk=0; kk<2*(sbporder-1)+1; kk++) {
                for (int kkk=0; kkk<ndim; kkk++) {
                    fdiff += (fdcoeff[0][kk]*jac[i][j][k-sbporder+1+kk]*
                              m[kkk][i][j][k-sbporder+1+kk]*
                              f[findex[kkk]][i][j][k-sbporder+1+kk]);
                }
            }
        }
        break;
    default:
        cerr << "Error in calculating finite difference, invalid direction " << dir << "\n";
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    
    return fdiff;
    
}


double fd_type::nonc(double*** f, const int i, const int j, const int k, const coord c, const int dir) const {
    // returns non-conservative finite difference of the given field for the indices given
    // in the specified direction
    
    assert(dir == 0 || dir == 1 || dir == 2);
    assert(i >=0 && i < c.get_nx_tot(0));
    assert(j >=0 && j < c.get_nx_tot(1));
    assert(k >=0 && k < c.get_nx_tot(2));
    
    double fdiff = 0.;
	int n, n_loc, m_loc, m_ghost;
    
    if (dir == 0) {
		n = c.get_nx(0);
		n_loc = c.get_nx_loc(0);
		m_loc = c.get_xm_loc(0)-c.get_xm(0);
        m_ghost = c.get_xm_ghost(0);
        if (i<2*(sbporder-1) && m_loc == 0) {
            for (int ii=0; ii<3*(sbporder-1); ii++) {
                fdiff += fdcoeff[i+1][ii]*f[ii][j][k];
            }
        } else if (i>n_loc+m_ghost-1-2*(sbporder-1) && m_loc+n_loc == n) {
            for (int ii=0; ii<3*(sbporder-1); ii++) {
                fdiff -= fdcoeff[n_loc+m_ghost-i][ii]*f[n_loc+m_ghost-1-ii][j][k];
            }
        } else {
            for (int ii=0; ii<2*(sbporder-1)+1; ii++) {
                fdiff += fdcoeff[0][ii]*f[i-sbporder+1+ii][j][k];
            }
        }
    } else if (dir == 1) {
		n = c.get_nx(1);
		n_loc = c.get_nx_loc(1);
		m_loc = c.get_xm_loc(1)-c.get_xm(1);
        m_ghost = c.get_xm_ghost(1);
        if (j<2*(sbporder-1) && m_loc == 0) {
            for (int jj=0; jj<3*(sbporder-1); jj++) {
                fdiff += fdcoeff[j+1][jj]*f[i][jj][k];
            }
        } else if (j>n_loc+m_ghost-1-2*(sbporder-1) && m_loc+n_loc == n) {
            for (int jj=0; jj<3*(sbporder-1); jj++) {
                fdiff -= fdcoeff[n_loc+m_ghost-j][jj]*f[i][n_loc+m_ghost-1-jj][k];
            }
        } else {
            for (int jj=0; jj<2*(sbporder-1)+1; jj++) {
                fdiff += fdcoeff[0][jj]*f[i][j-sbporder+1+jj][k];
            }
        }
    } else if (dir == 2) {
		n = c.get_nx(2);
		n_loc = c.get_nx_loc(2);
		m_loc = c.get_xm_loc(2)-c.get_xm(2);
        m_ghost = c.get_xm_ghost(2);
        if (k<2*(sbporder-1) && m_loc == 0) {
            for (int kk=0; kk<3*(sbporder-1); kk++) {
                fdiff += fdcoeff[k+1][kk]*f[i][j][kk];
            }
        } else if (k>n_loc+m_ghost-1-2*(sbporder-1) && m_loc+n_loc == n) {
            for (int kk=0; kk<3*(sbporder-1); kk++) {
                fdiff -= fdcoeff[n_loc+m_ghost-k][kk]*f[i][j][n_loc+m_ghost-1-kk];
            }
        } else {
            for (int kk=0; kk<2*(sbporder-1)+1; kk++) {
                fdiff += fdcoeff[0][kk]*f[i][j][k-sbporder+1+kk];
            }
        }
    } else {
        cerr << "Error in calculating finite difference, invalid direction " << dir << "\n";
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    
    return fdiff;

}

double fd_type::diss(double*** f, const int i, const int j, const int k, const coord c, const int dir) const {
    // returns disscoeff given indices
    
    assert(dir == 0 || dir == 1 || dir == 2);
    assert(i >=0 && i < c.get_nx_tot(0));
    assert(j >=0 && j < c.get_nx_tot(1));
    assert(k >=0 && k < c.get_nx_tot(2));
    
    double diss = 0.;
	int n, n_loc, m_loc, m_ghost;
    
    if (dir == 0) {
		n = c.get_nx(0);
		n_loc = c.get_nx_loc(0);
		m_loc = c.get_xm_loc(0);
        m_ghost = c.get_xm_ghost(0);
        if (i<2*(sbporder-1) && m_loc == 0) {
            for (int ii=0; ii<3*(sbporder-1); ii++) {
                diss += fdcoeff[i+1][ii]*f[ii][j][k];
            }
        } else if (i>n_loc+m_ghost-1-2*(sbporder-1) && m_loc+n_loc == n) {
            for (int ii=0; ii<3*(sbporder-1); ii++) {
                diss += fdcoeff[n_loc+m_ghost-i][ii]*f[n_loc+m_ghost-1-ii][j][k];
            }
        } else {
            for (int ii=0; ii<2*(sbporder-1)+1; ii++) {
                diss += fdcoeff[0][ii]*f[i-sbporder+1+ii][j][k];
            }
        }
    } else if (dir == 1) {
		n = c.get_nx(1);
		n_loc = c.get_nx_loc(1);
		m_loc = c.get_xm_loc(1);
        m_ghost = c.get_xm_ghost(1);
        if (j<2*(sbporder-1) && m_loc == 0) {
            for (int jj=0; jj<3*(sbporder-1); jj++) {
                diss += fdcoeff[j+1][jj]*f[i][jj][k];
            }
        } else if (j>n_loc-1-2*(sbporder-1) && m_loc+n_loc == n) {
            for (int jj=0; jj<3*(sbporder-1); jj++) {
                diss += fdcoeff[n_loc+m_ghost-j][jj]*f[i][n_loc+m_ghost-1-jj][k];
            }
        } else {
            for (int jj=0; jj<2*(sbporder-1)+1; jj++) {
                diss += fdcoeff[0][jj]*f[i][j-sbporder+1+jj][k];
            }
        }
    } else if (dir == 2) {
		n = c.get_nx(2);
		n_loc = c.get_nx_loc(2);
		m_loc = c.get_xm_loc(2);
        m_ghost = c.get_xm_ghost(2);
        if (k<2*(sbporder-1) && m_loc == 0) {
            for (int kk=0; kk<3*(sbporder-1); kk++) {
                diss += fdcoeff[k+1][kk]*f[i][j][kk];
            }
        } else if (k>n_loc+m_ghost-1-2*(sbporder-1) && m_loc+n_loc-1 == n) {
            for (int kk=0; kk<3*(sbporder-1); kk++) {
                diss += fdcoeff[n_loc+m_ghost-k][kk]*f[i][j][n_loc+m_ghost-1-kk];
            }
        } else {
            for (int kk=0; kk<2*(sbporder-1)+1; kk++) {
                diss += fdcoeff[0][kk]*f[i][j][k-sbporder+1+kk];
            }
        }
    } else {
        cerr << "Error in calculating dissipation, invalid direction " << dir << "\n";
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    
    return diss;
}