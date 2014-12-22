/************************************************************
  $Header: /home/dieter/sweph/RCS/swephexp.h,v 1.75 2009/04/08 07:19:08 dieter Exp $
  SWISSEPH: exported definitions and constants 

  This file represents the standard application interface (API)
  to the Swiss Ephemeris.

  A C programmer needs only to include this file, and link his code
  with the SwissEph library.

  The function calls are documented in the Programmer's documentation,
  which is online in HTML format.

  Structure of this file:
    Public API definitions
    Internal developer's definitions
    Public API functions.

  Authors: Dieter Koch and Alois Treindl, Astrodienst Zurich

************************************************************/
/* Copyright (C) 1997 - 2008 Astrodienst AG, Switzerland.  All rights reserved.

  License conditions
  ------------------

  This file is part of Swiss Ephemeris.

  Swiss Ephemeris is distributed with NO WARRANTY OF ANY KIND.  No author
  or distributor accepts any responsibility for the consequences of using it,
  or for whether it serves any particular purpose or works at all, unless he
  or she says so in writing.  

  Swiss Ephemeris is made available by its authors under a dual licensing
  system. The software developer, who uses any part of Swiss Ephemeris
  in his or her software, must choose between one of the two license models,
  which are
  a) GNU public license version 2 or later
  b) Swiss Ephemeris Professional License

  The choice must be made before the software developer distributes software
  containing parts of Swiss Ephemeris to others, and before any public
  service using the developed software is activated.

  If the developer choses the GNU GPL software license, he or she must fulfill
  the conditions of that license, which includes the obligation to place his
  or her whole software project under the GNU GPL or a compatible license.
  See http://www.gnu.org/licenses/old-licenses/gpl-2.0.html

  If the developer choses the Swiss Ephemeris Professional license,
  he must follow the instructions as found in http://www.astro.com/swisseph/ 
  and purchase the Swiss Ephemeris Professional Edition from Astrodienst
  and sign the corresponding license contract.

  The License grants you the right to use, copy, modify and redistribute
  Swiss Ephemeris, but only under certain conditions described in the License.
  Among other things, the License requires that the copyright notices and
  this notice be preserved on all copies.

  Authors of the Swiss Ephemeris: Dieter Koch and Alois Treindl

  The authors of Swiss Ephemeris have no control or influence over any of
  the derived works, i.e. over software or services created by other
  programmers which use Swiss Ephemeris functions.

  The names of the authors or of the copyright holder (Astrodienst) must not
  be used for promoting any software, product or service which uses or contains
  the Swiss Ephemeris. This copyright notice is the ONLY place where the
  names of the authors can legally appear, except in cases where they have
  given special permission in writing.

  The trademarks 'Swiss Ephemeris' and 'Swiss Ephemeris inside' may be used
  for promoting such software, products or services.
*/
import std.string;
import std.math;
import std.stdio;
import std.conv;


enum  MY_TRUE= 1;  /* for use in other defines, before TRUE is defined */
enum  MY_FALSE=0;  /* for use in other defines, before TRUE is defined */


enum INTEL_BYTE_ORDER=1;
alias int32=int;
alias in64=long;
alias int16=short;
alias uint32=uint;
alias REAL8=double;
alias int4=int;
alias uint4=int;
alias AS_BOOL=int;
alias UINT2=short;

enum TRUE=1;
enum FALSE=0;
enum OK=0;
enum ERR=-1;
alias UCHAR=ubyte;
/*
#define UCP (UCHAR*)
#define SCP (char*)
*/
immutable ODEGREE_STRING= "°";  /* degree as string, utf8 encoding */
 
immutable HUGE= 1.7E+308;     /* biggest value for REAL8 */
immutable M_PI= 3.14159265358979323846;
 
//#define forward static

enum AS_MAXCH= 256;    /* used for string declarations, allowing 255 char+\0 */
 
immutable DEGTORAD=0.0174532925199433;
immutable RADTODEG=57.2957795130823;
 
alias centisec=int32;
alias CSEC=centisec;

immutable DEG=360000;  /* degree expressed in centiseconds */
immutable DEG7_30=(2700000);  /* 7.5 degrees */
immutable DEG15   =(15 * DEG);
immutable DEG24   =(24 * DEG);
immutable DEG30   =(30 * DEG);
immutable DEG60   =(60 * DEG);
immutable DEG90   =(90 * DEG);
immutable DEG120  =(120 * DEG);
immutable DEG150  =(150 * DEG);
immutable DEG180  =(180 * DEG);
immutable DEG270  =(270 * DEG);
immutable DEG360  =(360 * DEG);
 
immutable CSTORAD =4.84813681109536E-08; /* centisec to rad: pi / 180 /3600/100 */
immutable RADTOCS  =2.06264806247096E+07; /* rad to centisec 180*3600*100/pi */
 
immutable CS2DEG=(1.0/360000.0);  /* centisec to degree */
immutable BFILE_R_ACCESS="r"; /* open binary file for reading */
immutable BFILE_RW_ACCESS="r+"; /* open binary file for writing and reading */
immutable BFILE_W_CREATE="w"; /* create/open binary file for write*/
immutable BFILE_A_ACCESS="a+";  /* create/open binary file for append*/
immutable FILE_R_ACCESS="r";  /* open text file for reading */
immutable FILE_RW_ACCESS="r+";  /* open text file for writing and reading */
immutable FILE_W_CREATE="w";  /* create/open text file for write*/
immutable FILE_A_ACCESS="a+"; /* create/open text file for append*/
immutable O_BINARY=0;   /* for open(), not defined in Unix */
immutable OPEN_MODE=std.conv.octal!666; /* default file creation mode */
immutable DIR_GLUE="/";   /* glue string for directory/file */
immutable PATH_SEPARATOR=";:";  /* semicolon or colon may be used */


/*
#    define FAR far
#    define MALLOC _fmalloc
#    define CALLOC _fcalloc
#    define FREE _ffree */

enum eph_flg
{
   jpleph =1, /* use JPL ephemeris */
   swisseph =2, /* use SWISSEPH ephemeris */
   mosheph =4, /* use Moshier ephemeris */
   helio =8, /* return heliocentric position */
   TRUEPOS =16, /* return true positions, not apparent */
   J2000 =32, /* no precession, i.e. give J2000 equinox */
   NONUT =64, /* no nutation, i.e. mean equinox of date */
   SPEED3 =128, /* speed from 3 positions (do not use it, */
   SPEED =256, /* high precision speed */ 
   NOGDEFL =512, /* turn off gravitational deflection */
   NOABERR =1024, /* turn off 'annual' aberration of light */
   EQUATORIAL =(2*1024), /* equatorial positions are wanted */
   XYZ =(4*1024), /* cartesian, not polar, coordinates */
   RADIANS =(8*1024), /* coordinates in radians, not degrees */
   bary =(16*1024), /* barycentric positions */
   topo =(32*1024), /* topocentric positions */
   sidereal =(64*1024), /* sidereal positions */
   ICRS =(128*1024) ,/* ICRS (DE406 reference frame) */
   DPSIDEPS_1980 =(256*1024), /* reproduce JPL Horizons */
   JPLHOR =DPSIDEPS_1980,
   JPLHOR_APPROX =(512*1024), /* approximate JPL Horizons 1962 - today */
   DEFAULTEPH =swisseph,
}

void throwOnError(int x)
{

}

struct STAR
{
  string name;
  double longitude;
  double latitude;
  double distance;
  double speed_long;
  double speed_lat;
  double speed_dist;
}

struct HOUSES
{
  double[13] cusps;
  double[10] ascmc;
}

struct EPHEMERIS
{
  char[255] serr;
  string path;

  enum SE_JUL_CAL  =0;
  enum GREG_CAL =1;
  enum ECL_NUT     =-1;      
  enum SUN =0;
  enum MOON =1;
  enum MERCURY =2;
  enum VENUS =3;
  enum MARS =4;
  enum JUPITER =5;
  enum SATURN =6;
  enum URANUS =7;
  enum NEPTUNE =8;
  enum PLUTO =9;
  enum MEAN_NODE =10;
  enum TRUE_NODE =11;
  enum MEAN_APOG =12;
  enum OSCU_APOG =13;
  enum EARTH =14;
  enum CHIRON =15;
  enum PHOLUS =16;
  enum CERES =17;
  enum PALLAS =18;
  enum JUNO =19;
  enum VESTA =20;
  enum INTP_APOG =21;
  enum INTP_PERG =22;
  enum NPLANETS =23;
  enum AST_OFFSET =10000;
  enum VARUNA =(AST_OFFSET+20000);
  enum FICT_OFFSET =  40;
  enum FICT_OFFSET_1 =  39;
  enum FICT_MAX = 999;
  enum NFICT_ELEM =15;
  enum COMET_OFFSET =1000;
  enum NALL_NAT_POINTS =(NPLANETS+NFICT_ELEM);

  /* Hamburger or Uranian "planets" */
  enum CUPIDO = 40;
  enum HADES =  41;
  enum ZEUS = 42;
  enum KRONOS = 43;
  enum APOLLON =  44;
  enum ADMETOS =  45;
  enum VULKANUS = 46;
  enum POSEIDON = 47;
  /* other fictitious bodies */
  enum ISIS = 48;
  enum NIBIRU = 49;
  enum HARRINGTON =50;
  enum NEPTUNE_LEVERRIER =51;
  enum NEPTUNE_ADAMS =52;
  enum PLUTO_LOWELL =53;
  enum PLUTO_PICKERING =54;
  enum VULCAN =   55;
  enum WHITE_MOON =   56;
  enum PROSERPINA =   57;
  enum WALDEMATH =    58;
  enum FIXSTAR =-10;

  enum ASC      =0;
  enum MC     =1;
  enum ARMC     =2;
  enum VERTEX   =3;
  enum EQUASC =   4; /*"equatorialascendant"*/;
  enum COASC1   =5; /* ="co-ascendant"(W.Koch)*/;
  enum COASC2   =6; /* ="co-ascendant"(M.Munkasey)*/;
  enum POLASC   =7; /* ="polarascendant"(M.Munkasey)*/;
  enum NASCMC   =8;

  /*
   * flag bits for parameter iflag in function swe_calc()
   * The flag bits are defined in such a way that iflag = 0 delivers what one
   * usually wants:
   *    - the default ephemeris (SWISS EPHEMERIS) is used,
   *    - apparent geocentric positions referring to the true equinox of date
   *      are returned.
   * If not only coordinates, but also speed values are required, use 
   * flag = SEFLG_SPEED.
   *
   * The 'L' behind the number indicates that 32-bit integers (Long) are used.
   */

  enum SIDBITS =256 ;
  /* for projection onto ecliptic of t0 */
  enum SIDBIT_ECL_T0 =256 ;
  /* for projection onto solar system plane */
  enum SIDBIT_SSY_PLANE =512 ;

  /* sidereal modes (ayanamsas) */
  enum SIDM_FAGAN_BRADLEY =0 ;
  enum SIDM_LAHIRI =1 ;
  enum SIDM_DELUCE =2 ;
  enum SIDM_RAMAN =3 ;
  enum SIDM_USHASHASHI =4 ;
  enum SIDM_KRISHNAMURTI =5 ;
  enum SIDM_DJWHAL_KHUL =6 ;
  enum SIDM_YUKTESHWAR =7 ;
  enum SIDM_JN_BHASIN =8 ;
  enum SIDM_BABYL_KUGLER1 =9 ;
  enum SIDM_BABYL_KUGLER2 =10 ;
  enum SIDM_BABYL_KUGLER3 =11 ;
  enum SIDM_BABYL_HUBER =12 ;
  enum SIDM_BABYL_ETPSC =13 ;
  enum SIDM_ALDEBARAN_15TAU =14 ;
  enum SIDM_HIPPARCHOS =15 ;
  enum SIDM_SASSANIAN =16 ;
  enum SIDM_GALCENT_0SAG =17 ;
  enum SIDM_J2000 =18 ;
  enum SIDM_J1900 =19 ;
  enum SIDM_B1950 =20 ;
  enum SIDM_SURYASIDDHANTA =21 ;
  enum SIDM_SURYASIDDHANTA_MSUN =22 ;
  enum SIDM_ARYABHATA =23 ;
  enum SIDM_ARYABHATA_MSUN =24 ;
  enum SIDM_SS_REVATI =25 ;
  enum SIDM_SS_CITRA =26 ;
  enum SIDM_TRUE_CITRA =27 ;
  enum SIDM_TRUE_REVATI =28 ;
  enum SIDM_USER =255 ;

  enum NSIDM_PREDEF =29 ;

  /* used for swe_nod_aps(): */
  enum NODBIT_MEAN =1 /* mean nodes/apsides */ ;
  enum NODBIT_OSCU =2 /* osculating nodes/apsides */ ;
  enum NODBIT_OSCU_BAR =4 /* same, but motion about solar system barycenter is considered */ ;
  enum NODBIT_FOPOINT =256 /* focal point of orbit instead of aphelion */ ;


  enum MAX_STNAME =256; /* maximum size of fixstar name; ;
                                           * the parameter star in swe_fixstar
             * must allow twice this space for
                   * the returned star name.
             */

  /* defines for eclipse computations */

  enum ECL_CENTRAL =1 ;
  enum ECL_NONCENTRAL =2 ;
  enum ECL_TOTAL =4 ;
  enum ECL_ANNULAR =8 ;
  enum ECL_PARTIAL =16 ;
  enum ECL_ANNULAR_TOTAL =32 ;
  enum ECL_PENUMBRAL =64 ;
  enum ECL_ALLTYPES_SOLAR =(ECL_CENTRAL|ECL_NONCENTRAL|ECL_TOTAL|ECL_ANNULAR|ECL_PARTIAL|ECL_ANNULAR_TOTAL) ;
  enum ECL_ALLTYPES_LUNAR =(ECL_TOTAL|ECL_PARTIAL|ECL_PENUMBRAL) ;
  enum ECL_VISIBLE =128 ;
  enum ECL_MAX_VISIBLE =256 ;
  enum ECL_1ST_VISIBLE =512 /* begin of partial eclipse */ ;
  enum ECL_PARTBEG_VISIBLE =512 /* begin of partial eclipse */ ;
  enum ECL_2ND_VISIBLE =1024 /* begin of total eclipse */ ;
  enum ECL_TOTBEG_VISIBLE =1024 /* begin of total eclipse */ ;
  enum ECL_3RD_VISIBLE =2048 /* end of total eclipse */ ;
  enum ECL_TOTEND_VISIBLE =2048 /* end of total eclipse */ ;
  enum ECL_4TH_VISIBLE =4096 /* end of partial eclipse */ ;
  enum ECL_PARTEND_VISIBLE =4096 /* end of partial eclipse */ ;
  enum ECL_PENUMBBEG_VISIBLE =8192 /* begin of penumbral eclipse */ ;
  enum ECL_PENUMBEND_VISIBLE =16384 /* end of penumbral eclipse */ ;
  enum ECL_OCC_BEG_DAYLIGHT =8192 /* occultation begins during the day */ ;
  enum ECL_OCC_END_DAYLIGHT =16384 /* occultation ends during the day */ ;
  enum ECL_ONE_TRY =(32*1024) ;

      /* check if the next conjunction of the moon with
       * a planet is an occultation; don't search further */

  /* for swe_ritransit() */
  enum CALC_RISE =1 ;
  enum CALC_SET =2 ;
  enum CALC_MTRANSIT =4 ;
  enum CALC_ITRANSIT =8 ;
  enum BIT_DISC_CENTER =256 /* to be or'ed to CALC_RISE/SET */ ;
  enum BIT_DISC_BOTTOM =8192 /* to be or'ed to CALC_RISE/SET*/ ;
  enum BIT_NO_REFRACTION =512 /* to be or'ed to CALC_RISE/SET*/ ;
  enum BIT_CIVIL_TWILIGHT =1024 /* to be or'ed to CALC_RISE/SET */ ;
  enum BIT_NAUTIC_TWILIGHT =2048 /* to be or'ed to CALC_RISE/SET */ ;
  enum BIT_ASTRO_TWILIGHT =4096 /* to be or'ed to CALC_RISE/SET */ ;
  enum BIT_FIXED_DISC_SIZE =(16*1024) /* or'ed to CALC_RISE/SET: */ ;
  enum ECL2HOR =0 ;
  enum EQU2HOR =1 ;
  enum HOR2ECL =0 ;
  enum HOR2EQU =1 ;
  enum TRUE_TO_APP =0 ;
  enum APP_TO_TRUE =1 ;
  enum DE_NUMBER =431 ;

  enum FNAME_DE200 ="de200.eph" ;
  enum FNAME_DE403 ="de403.eph" ;
  enum FNAME_DE404 ="de404.eph" ;
  enum FNAME_DE405 ="de405.eph" ;
  enum FNAME_DE406 ="de406.eph" ;
  enum FNAME_DE431 ="de431.eph" ;
  enum FNAME_DFT =FNAME_DE431 ;
  enum FNAME_DFT2 =FNAME_DE406 ;
  enum STARFILE_OLD ="fixstars.cat" ;
  enum STARFILE ="sefstars.txt" ;
  enum ASTNAMFILE ="seasnam.txt" ;
  enum FICTFILE ="seorbel.txt" ;
  enum HELIACAL_RISING =1 ;
  enum HELIACAL_SETTING =2 ;
  enum MORNING_FIRST =HELIACAL_RISING ;
  enum EVENING_LAST =HELIACAL_SETTING ;
  enum EVENING_FIRST =3 ;
  enum MORNING_LAST =4 ;
  enum ACRONYCHAL_RISING =5 /* still not implemented */ ;
  enum ACRONYCHAL_SETTING =6 /* still not implemented */ ;
  enum COSMICAL_SETTING =ACRONYCHAL_SETTING ;
  enum HELFLAG_LONG_SEARCH =128 ;
  enum HELFLAG_HIGH_PRECISION =256 ;
  enum HELFLAG_OPTICAL_PARAMS =512 ;
  enum HELFLAG_NO_DETAILS =1024 ;
  enum HELFLAG_SEARCH_1_PERIOD =(1 << 11) /* 2048 */ ;
  enum HELFLAG_VISLIM_DARK =(1 << 12) /* 4096 */ ;
  enum HELFLAG_VISLIM_NOMOON =(1 << 13) /* 8192 */ ;
  enum HELFLAG_VISLIM_PHOTOPIC =(1 << 14) /* 16384 */ ;
  enum HELFLAG_AV =(1 << 15) /* 32768 */ ;
  enum HELFLAG_AVKIND_VR =(1 << 15) /* 32768 */ ;
  enum HELFLAG_AVKIND_PTO =(1 << 16) ;
  enum HELFLAG_AVKIND_MIN7 =(1 << 17) ;
  enum HELFLAG_AVKIND_MIN9 =(1 << 18) ;


  enum EPHE_PATH = ".:/users/ephe2/:/users/ephe/";
  enum SPLIT_DEG_ROUND_SEC   = 1;
  enum SPLIT_DEG_ROUND_MIN   = 2;
  enum SPLIT_DEG_ROUND_DEG   = 4;
  enum SPLIT_DEG_ZODIACAL    = 8;
  enum SPLIT_DEG_KEEP_SIGN   =16; /* don't round to next sign, 
             * e.g. 29.9999999 will be rounded
             * to 29d59'59" (or 29d59' or 29d) */
  enum SPLIT_DEG_KEEP_DEG    =32; /* don't round to next degree
             * e.g. 13.9999999 will be rounded
             * to 13d59'59" (or 13d59' or 13d) */


  enum HELFLAG_AVKIND =(HELFLAG_AVKIND_VR|HELFLAG_AVKIND_PTO|HELFLAG_AVKIND_MIN7|HELFLAG_AVKIND_MIN9) ;
  enum TJD_INVALID =99999999.0 ;
  enum SIMULATE_VICTORVB =1 ;

  enum HELIACAL_LONG_SEARCH =128 ;
  enum HELIACAL_HIGH_PRECISION =256 ;
  enum HELIACAL_OPTICAL_PARAMS =512 ;
  enum HELIACAL_NO_DETAILS =1024 ;
  enum HELIACAL_SEARCH_1_PERIOD =(1 << 11) /* 2048 */ ;
  enum HELIACAL_VISLIM_DARK =(1 << 12) /* 4096 */ ;
  enum HELIACAL_VISLIM_NOMOON =(1 << 13) /* 8192 */ ;
  enum HELIACAL_VISLIM_PHOTOPIC =(1 << 14) /* 16384 */ ;
  enum HELIACAL_AVKIND_VR =(1 << 15) /* 32768 */ ;
  enum HELIACAL_AVKIND_PTO =(1 << 16) ;
  enum HELIACAL_AVKIND_MIN7 =(1 << 17) ;
  enum HELIACAL_AVKIND_MIN9 =(1 << 18) ;
  enum HELIACAL_AVKIND =(HELFLAG_AVKIND_VR|HELFLAG_AVKIND_PTO|HELFLAG_AVKIND_MIN7|HELFLAG_AVKIND_MIN9) ;
  enum PHOTOPIC_FLAG =0 ;
  enum SCOTOPIC_FLAG =1 ;
  enum MIXEDOPIC_FLAG =2 ;

  string ZtoString(char[] c)
  {
    return to!string(fromStringz(cast(char*)c));
  }

  string ZtoString(char* c)
  {
    return to!string(fromStringz(c));
  }

  void set_ephe_path(string path)
  {
    swe_set_ephe_path(toStringz(path));
  }
  
  void set_jpl_file(string fname)
  {
    swe_set_jpl_file(toStringz(fname));
  }

  int heliacal_ut(double tjdstart_ut, double *geopos, double *datm, double *dobs, string ObjectName, int TypeEvent, int32 iflag, double *dret, ref string serr)
  {
    int ret;
    char[1024] buf;
    ret=swe_heliacal_ut(tjdstart_ut, geopos, datm, dobs, toStringz(ObjectName), TypeEvent, iflag, dret, cast(char*)buf);
    serr=to!string(buf);
    return ret;
  }

  int heliacal_pheno_ut(double tjd_ut, double *geopos, double *datm, double *dobs, string ObjectName, int TypeEvent, int helflag, double *darr, ref string serr)
  {
    int ret;
    char[1024] buf;
    ret=swe_heliacal_pheno_ut(tjd_ut, geopos, datm, dobs, toStringz(ObjectName),  TypeEvent,  helflag, darr, cast(char*)buf);
    serr=to!string(buf);
    return ret; 
  }

  int vis_limit_mag(double tjdut, double *geopos, double *datm, double *dobs, string ObjectName, int helflag, double *dret,ref string serr)
  {
    int ret;
    char[1024] buf;
    ret=swe_vis_limit_mag(tjdut, geopos,datm, dobs, toStringz(ObjectName), helflag, dret,cast(char*)buf);
    serr=to!string(buf);
    return ret; 
  }

  int heliacal_angle(double tjdut, double *dgeo, double *datm, double *dobs, int helflag, double mag, double azi_obj, double azi_sun, double azi_moon, double alt_moon, double *dret, ref string serr)
  {
    int ret;
    char[1024] buf;
    ret=swe_heliacal_angle(tjdut,dgeo, datm, dobs, helflag, mag, azi_obj, azi_sun,azi_moon,  alt_moon, dret, cast(char*)buf);
    serr=to!string(buf);
    return ret; 
  }

  int topo_arcus_visionis(double tjdut, double *dgeo, double *datm, double *dobs, int32 helflag, double mag, double azi_obj, double alt_obj, double azi_sun, double azi_moon, double alt_moon, double *dret, ref string serr)
  {
    int ret;
    char[1024] buf;
    ret=swe_topo_arcus_visionis(tjdut,dgeo,datm,dobs,helflag, mag,azi_obj,alt_obj,azi_sun,azi_moon,alt_moon, dret, cast(char*)buf);
    serr=to!string(buf);
    return ret; 
  }

  string get_version()
  {
    char[1024] buf;
    char* ptr;
    ptr=swe_version(cast(char*)buf);
    return to!string(ptr);
  }

  string get_version(ref string verbuf)
  {
    char[1024] buf=verbuf;
    char* ptr;
    ptr=swe_version(cast(char*)buf);
    verbuf=to!string(buf);
    return verbuf;
  }

  int calc(double tjd, int ipl, int iflag, ref double[6] xx, ref string serr)
  {
    int ret;
    char[255] buf;
    ret=swe_calc(tjd, ipl, iflag, cast(double*)xx, cast(char*)buf);
    serr=to!string(buf);
    return ret; 
  }

  double[6] calc_ut(double tjd_ut, int ipl, int iflag)
  {
    int ret;
    double[6] xx;
    throwOnError(swe_calc_ut(tjd_ut, ipl, iflag, cast(double*)xx, cast(char*)serr));
    return xx; 
  }

  string strerr()
  {
    return ZtoString(serr);
  }

  
  STAR fixstar(string star, double tjd, int iflag)
  {
    STAR ret;
    double[6] xx;
    char[1024] buf,starbuf;
    starbuf=star;
    throwOnError(swe_fixstar(cast(char*)starbuf, tjd,iflag, cast(double*)xx, cast(char *)serr));
    ret.name=ZtoString(starbuf);
    ret.longitude=xx[0];
    ret.latitude=xx[1];
    ret.distance=xx[2];
    ret.speed_long=xx[3];
    ret.speed_lat=xx[4];
    ret.speed_dist=xx[5];
    return ret; 
  }


  STAR fixstar_ut(string star, double tjd_ut, int32 iflag)
  {
    STAR ret;
    double[6] xx;
    char[255] starbuf;
    starbuf=star;
    throwOnError(swe_fixstar_ut(cast(char*)starbuf,  tjd_ut, iflag, cast(double*)xx, cast(char *)serr));
    ret.name=ZtoString(starbuf);
    ret.longitude=xx[0];
    ret.latitude=xx[1];
    ret.distance=xx[2];
    ret.speed_long=xx[3];
    ret.speed_lat=xx[4];
    ret.speed_dist=xx[5];
    return ret;
  }

  double fixstar_mag(string star)
  {
    double mag;
    int ret;
    char[1024] starbuf;
    starbuf=star;
    throwOnError(swe_fixstar_mag(cast(char*)starbuf,  &mag,cast(char*)serr));
    return mag;
  }

  string get_planet_name(int ipl)
  {
    int ret;
    char[255] buf;
    buf[0]=0;
    auto retname=swe_get_planet_name(ipl,cast(char*)buf);
    return to!string(fromStringz(cast(char*)buf));
  }

  string get_ayanamsa_name(int isidmode)
  {
    return to!string(fromStringz(swe_get_ayanamsa_name(isidmode)));
  }

  int utc_to_jd(int32 iyear, int32 imonth, int32 iday, int32 ihour, int32 imin, double dsec, int32 gregflag,double *dret, ref string serr)
  {
    int ret;
    char[1024] buf;
    ret= swe_utc_to_jd(iyear, imonth, iday, ihour, imin, dsec, gregflag, dret,cast(char*) serr);
    serr=to!string(buf);
    return ret;
  }

  HOUSES houses(double tjd_ut, double geolat, double geolon, int hsys)
  {
    HOUSES ret;
    throwOnError(swe_houses(tjd_ut, geolat, geolon, hsys, cast(double*)ret.cusps, cast(double*)ret.ascmc));
    return ret;
  }

  HOUSES houses_ex(double tjd_ut, int32 iflag, double geolat, double geolon, int hsys)
  {
    HOUSES ret;
    throwOnError(swe_houses_ex(tjd_ut,iflag,geolat,geolon,hsys,cast(double*)ret.cusps,cast(double*)ret.ascmc));
    return ret;
  }

  HOUSES houses_armc(double armc, double geolat, double eps, int hsys)
  {
    HOUSES ret;
    throwOnError(swe_houses_armc(armc, geolat, eps, hsys,cast(double*)ret.cusps,cast(double*)ret.ascmc));
    return ret;
  }

  double house_pos(double armc, double geolat, double eps, int hsys, in double[2] xpin)
  {
    double ret;
    return swe_house_pos(armc,geolat,eps,hsys,cast(double*)xpin, cast(char*)serr);
  }

  double house_pos(double armc, double geolat, double eps, int hsys, double ecl_longitude, double ecl_latitude)
  {
    double ret;
    double[2] xpin;
    xpin[0]=ecl_longitude;
    xpin[1]=ecl_latitude;
    ret= swe_house_pos(armc,geolat,eps,hsys,cast(double*)xpin, cast(char*)serr);
    return ret;
  }

  string house_name(int hsys)
  {
    return to!string(swe_house_name(hsys));
  }

  int gauquelin_sector(double t_ut, int32 ipl, string starname, int iflag, int imeth, in double[3] geopos, double atpress, double attemp, out double dgsect, ref string serr)
  {
    int ret;
    char[1024] buf;
    ret=swe_gauquelin_sector(t_ut,ipl,toStringz(starname),iflag,imeth,cast(double*)geopos,atpress, attemp,cast(double*)&dgsect, cast(char*)buf);
    serr=to!string(buf);
    return ret;
  }

  int pheno(double tjd, int32 ipl, int32 iflag, out double[20] attr, ref string serr)
  {
    int ret;
    char[1024] buf;
    ret=swe_pheno(tjd,ipl,iflag,cast(double*)attr,cast(char*)buf);
    serr=to!string(buf);
    return ret;
  }

  int pheno_ut(double tjd_ut, int32 ipl, int32 iflag, out double[20] attr, ref string serr)
  {
    int ret;
    char[1024] buf;
    ret=swe_pheno_ut(tjd_ut,ipl,iflag,cast(double*)attr,cast(char*)buf);
    serr=to!string(buf);
    return ret; 
  }

  int time_equ(double tjd, ref double te, ref string serr)
  {
    int ret;
    char[1024] buf;
    ret=swe_time_equ( tjd, cast(double*)&te, cast(char*)buf);
    serr=to!string(buf);
    return ret; 
  }

  double[3] geopos(double longitude, double latitude, double height)
  {
    return[longitude,latitude,height];
  }

  int sol_eclipse_where(double tjd, int32 ifl, out double[2] geopos, out double[20] attr, ref string serr)
  {
    int ret;
    char[1024] buf;
    ret=swe_sol_eclipse_where(tjd,ifl, cast(double*)geopos, cast(double*)attr, cast(char*) buf);
    serr=to!string(buf);
    return ret; 
  }

  int sol_eclipse_how(double tjd, int32 ifl, in double[3] geopos, out double[20] attr, ref string serr)
  {
    int ret;
    char[1024] buf;
    ret=swe_sol_eclipse_where(tjd,ifl, cast(double*)geopos, cast(double*)attr, cast(char*) buf);
    serr=to!string(buf);
    return ret; 
  }
  int sol_eclipse_when_loc(double tjd_start, int32 ifl, in double[3] geopos, out double[10] tret, out double[20] attr, bool backward, ref string serr)
  {
    int ret;
    char[1024] buf;
    int asbackward=(backward)?1:0;
    ret=swe_sol_eclipse_when_loc(tjd_start, ifl, cast(double*)geopos, cast(double*)tret, cast(double*)attr, asbackward, cast(char*)buf);
    serr=to!string(buf);
    return ret; 
  }

  int sol_eclipse_when_glob(double tjd_start, int32 ifl, int ifltype, out double[10] tret, bool backward, ref string serr)
  {
    int ret;
    char[1024] buf;
    int asbackward=(backward)?1:0;
    ret=swe_sol_eclipse_when_glob( tjd_start, ifl,  ifltype, cast(double*)tret, asbackward, cast(char*)buf);
    serr=to!string(buf);
    return ret; 
  }

  int lun_occult_where(double tjd, int ipl, string starname, int32 ifl, out double[2] geopos, out double[20] attr, ref string serr)
  {
    int ret;
    char[1024] buf;
    ret=swe_lun_occult_where(tjd,ipl, toStringz(starname), ifl, cast(double*) geopos, cast(double*) attr, cast(char*)buf);
    serr=to!string(buf);
    return ret;  
  }
  int lun_occult_when_loc(double tjd_start, int ipl, string starname, int32 ifl, in double[3] geopos, out double[10] tret, out double[20] attr, bool backward, ref string serr)
  {
    int ret;
    char[1024] buf;
    int asbackward=(backward)?1:0;
    ret=swe_lun_occult_when_loc(tjd_start, ipl, toStringz(starname), ifl, cast(double*)geopos, cast(double*)tret, cast(double*) attr, asbackward, cast(char*)buf);
    serr=to!string(buf);
    return ret; 
  }

  int lun_occult_when_glob(double tjd_start, int32 ipl, string starname, int ifl, int ifltype, out double[10] tret, bool backward, ref string serr)
  {
    int ret;
    char[1024] buf;
    int asbackward=(backward)?1:0;
    ret=swe_lun_occult_when_glob(tjd_start, ipl, toStringz(starname),  ifl,  ifltype,cast(double*) tret,  asbackward, cast(char*)buf);
    serr=to!string(buf);
    return ret; 
  }

  int lun_eclipse_how(double tjd_ut, int32 ifl, in double[3] geopos, out double[20] attr, ref string serr)
  {
    int ret;
    char[1024] buf;
    ret= swe_lun_eclipse_how( tjd_ut,  ifl, cast(double*) geopos,cast(double*)attr, cast(char*)buf);
    serr=to!string(buf);
    return ret; 
  }
  int lun_eclipse_when(double tjd_start, int ifl, int ifltype, out double[10] tret, bool backward, ref string serr)
  {
    int ret;
    char[1024] buf;
    int asbackward=(backward)?1:0;
    ret=swe_sol_eclipse_when_glob( tjd_start, ifl,  ifltype, cast(double*)tret, asbackward, cast(char*)buf);
    serr=to!string(buf);
    return ret; 
  }
  /**
  int lun_eclipse_when_loc(double tjd_start, int32 ifl, in double[3] geopos, out double[10] tret,out double[20] attr, bool backward, ref string serr)
  {
    int ret;
    char[1024] buf;
    int asbackward=(backward)?1:0;
    ret=swe_lun_eclipse_when_loc( tjd_start, ifl,cast(double*) geopos, cast(double*) tret,cast(double*) attr, asbackward, cast(char*)buf);
    serr=to!string(buf);
    return ret; 
  }
  */
  int rise_trans_true_hor(double tjd_ut, int32 ipl, string starname, int epheflag, int rsmi, in double[3] geopos, double atpress, double attemp, double horhgt, double *tret, ref string serr)
  {
    int ret;
    char[1024] buf;
    ret= swe_rise_trans_true_hor(tjd_ut, ipl, toStringz(starname), epheflag,  rsmi, cast(double*)geopos,  atpress,  attemp,  horhgt, cast(double *)tret, cast(char*)buf);
    serr=to!string(buf);
    return ret; 
  }
  int rise_trans(double tjd_ut, int32 ipl, string starname, int epheflag, int rsmi, in double[3] geopos, double atpress, double attemp, double *tret, ref string serr)
  {
    int ret;
    char[1024] buf;
    ret= swe_rise_trans( tjd_ut,  ipl, toStringz(starname),  epheflag,  rsmi,cast(double*)geopos, atpress, attemp, tret, cast(char*)buf);
    serr=to!string(buf);
    return ret; 
  }
  int nod_aps(double tjd_et, int32 ipl, int32 iflag, int32  method, out double[6] xnasc, out double[6] xndsc, out double[6] xperi, out double[6] xaphe, ref string serr)
  {
    int ret;
    char[1024] buf;
    ret= swe_nod_aps(tjd_et, ipl, iflag,  method,cast(double*) xnasc,cast(double*)xndsc,cast(double*)xperi,cast(double*) xaphe, cast(char*)buf);
    serr=to!string(buf);
    return ret; 
  }

  int nod_aps_ut(double tjd_ut, int32 ipl, int32 iflag, int32  method, out double[6] xnasc, out double[6] xndsc, out double[6] xperi, out double[6] xaphe, ref string serr)
  {
    int ret;
    char[1024] buf;
    ret=swe_nod_aps_ut(tjd_ut, ipl, iflag, method, cast(double*) xnasc, cast(double*)xndsc, cast(double*) xperi, cast(double*) xaphe, cast(char*)buf);
    serr=to!string(buf);
    return ret; 
  }

  /**
  int swe_lmt_to_lat(double tjd_lmt, double geolon, ref double tjd_lat, ref string serr)
  {
    int ret;
    char[1024] buf;
    ret=swe_lmt_to_lat(tjd_lmt, geolon,cast(double*)&tjd_lat, cast(char*)buf);
    serr=to!string(buf);
    return ret; 
  }

  int swe_lat_to_lmt(double tjd_lat, double geolon, ref double tjd_lmt, ref string serr)
  {
    int ret;
    char[1024] buf;
    ret=swe_lat_to_lmt(tjd_lat,geolon,cast(double*)&tjd_lmt, cast(char*)buf);
    serr=to!string(buf);
    return ret; 
  }
  */
  string cs2timestr(CSEC t, int sep, AS_BOOL suppressZero, ref string a)
  {
    char abuf[1024]=a;
    string s=to!string(swe_cs2timestr(t, sep, suppressZero, cast(char*)abuf));
    a=to!string(abuf);
    return s;
  }

  string cs2lonlatstr(CSEC t, char pchar, char mchar, ref string s)
  {
    char sbuf[1024]=s;
    string sret=to!string(swe_cs2lonlatstr(t, pchar, mchar, cast(char*)sbuf));
    s=to!string(sbuf);
    return sret;
  }

  string cs2degstr(CSEC t, ref string a)
  {
    char sbuf[1024]=a;
    string sret=to!string(swe_cs2degstr(t,cast(char*)sbuf));
    a=to!string(sbuf);
    return sret;
  }

  void azalt(double tjd_ut, int calc_flag, in double[3] geopos, double atpress, double attemp, in double[3] xin, out double[3] xaz)
  {
    swe_azalt( tjd_ut, calc_flag, cast(double*) geopos, atpress,  attemp, cast(double *)xin, cast(double *)xaz);
  }

  void azalt_rev(double tjd_ut, int calc_flag, in double[3] geopos, in double[3] xin, out double[3] xout)
  {
    swe_azalt_rev(tjd_ut, calc_flag, cast(double *)geopos, cast(double *)xin, cast(double *)xout);
  }

  void cotrans(in double[3] xpo, out double[3] xpn, double eps)
  {
    swe_cotrans(cast(double*)xpo,cast(double*)xpn,eps);
  }
  void cotrans_sp(in double[6] xpo, out double[6] xpn, double eps)
  {
    swe_cotrans_sp(cast(double*)xpo,cast(double*)xpn,eps);
  }

  double julday(int year, int month, int day)
  {
    return swe_julday(year,month,day,0.0,1);
  }

  double julday(int year, int month, int day, double hour)
  {
    return swe_julday(year,month,day,hour,1);
  }

  double  difdegn (double p1, double p2)
  {
    return swe_difdegn(p1,p2);
  }
  centisec difcs2n(centisec p1, centisec p2)
  {
    return swe_difcs2n(p1,p2);
  }
  double difdeg2n(double p1, double p2)
  {
    return swe_difdeg2n(p1,p2);
  }
  double  difrad2n(double p1, double p2)
  {
    return swe_difrad2n(p1,p2);
  }
  centisec csroundsec(centisec x)
  {
    return swe_csroundsec(x);
  }
  int  d2l(double x)
  {
    return swe_d2l(x);
  }
  int  day_of_week(double jd)
  {
    return swe_day_of_week(jd);
  }
}

extern(C) int swe_gauquelin_sector(double t_ut, int32 ipl, const(char)*starname, int32 iflag, int32 imeth, double *geopos, double atpress, double attemp, double *dgsect, char *serr);
extern(C) int  swe_heliacal_ut(double tjdstart_ut, double *geopos, double *datm, double *dobs, immutable(char) *ObjectName, int32 TypeEvent, int32 iflag, double *dret, char *serr);
extern(C) int  swe_heliacal_pheno_ut(double tjd_ut, double *geopos, double *datm, double *dobs, immutable(char) *ObjectName, int32 TypeEvent, int32 helflag, double *darr, char *serr);
extern(C) int  swe_vis_limit_mag(double tjdut, double *geopos, double *datm, double *dobs, const(char *)ObjectName, int32 helflag, double *dret, char *serr);
extern(C) int  swe_heliacal_angle(double tjdut, double *dgeo, double *datm, double *dobs, int32 helflag, double mag, double azi_obj, double azi_sun, double azi_moon, double alt_moon, double *dret, char *serr);
extern(C) int  swe_topo_arcus_visionis(double tjdut, double *dgeo, double *datm, double *dobs, int32 helflag, double mag, double azi_obj, double alt_obj, double azi_sun, double azi_moon, double alt_moon, double *dret, char *serr);
extern(C) char * swe_version(char *);
extern(C) int  swe_calc(double tjd, int ipl, int32 iflag, double *xx, char *serr);
extern(C) int  swe_calc_ut(double tjd_ut, int32 ipl, int32 iflag, double *xx, char *serr);
extern(C) int  swe_fixstar(char *star, double tjd, int32 iflag, double *xx, char *serr);
extern(C) int  swe_fixstar_ut(char *star, double tjd_ut, int32 iflag, double *xx, char *serr);
extern(C) int  swe_fixstar_mag(char *star, double *mag, char *serr);
extern(C) void  swe_close();
extern(C) void  swe_set_ephe_path(const(char *)path);
extern(C) void  swe_set_jpl_file(const(char *)fname);
extern(C) char* swe_get_planet_name(int ipl, char *spname);
extern(C) void swe_set_topo(double geolon, double geolat, double geoalt);
extern(C) void swe_set_sid_mode(int32 sid_mode, double t0, double ayan_t0);
extern(C) double swe_get_ayanamsa(double tjd_et);
extern(C) double swe_get_ayanamsa_ut(double tjd_ut);
extern(C) char* swe_get_ayanamsa_name(int32 isidmode);
extern(C) int  swe_date_conversion(int y , int m , int d , double utime, char c, double *tjd);
extern(C) double  swe_julday(int year, int month, int day, double hour, int gregflag);
extern(C) void  swe_revjul (double jd, int gregflag, int *jyear, int *jmon, int *jday, double *jut);
extern(C) int  swe_utc_to_jd(int32 iyear, int32 imonth, int32 iday, int32 ihour, int32 imin, double dsec, int32 gregflag, double *dret, char *serr);
extern(C) void swe_jdet_to_utc(double tjd_et, int32 gregflag, int32 *iyear, int32 *imonth, int32 *iday, int32 *ihour, int32 *imin, double *dsec);
extern(C) void swe_jdut1_to_utc(double tjd_ut, int32 gregflag, int32 *iyear, int32 *imonth, int32 *iday, int32 *ihour, int32 *imin, double *dsec);
extern(C) void swe_utc_time_zone(int32 iyear, int32 imonth, int32 iday, int32 ihour, int32 imin, double dsec, double d_timezone, int32 *iyear_out, int32 *imonth_out, int32 *iday_out, int32 *ihour_out, int32 *imin_out, double *dsec_out);
extern(C) int  swe_houses(double tjd_ut, double geolat, double geolon, int hsys, double *cusps, double *ascmc);
extern(C) int  swe_houses_ex(double tjd_ut, int32 iflag, double geolat, double geolon, int hsys, double *cusps, double *ascmc);
extern(C) int  swe_houses_armc(double armc, double geolat, double eps, int hsys, double *cusps, double *ascmc);
extern(C) double swe_house_pos(double armc, double geolat, double eps, int hsys, double *xpin, char *serr);
extern(C) char * swe_house_name(int hsys);
extern(C) int  swe_gauquelin_sector(double t_ut, int32 ipl, const(char) *starname, int32 iflag, int32 imeth, const(double) *geopos, double atpress, double attemp, double *dgsect, char *serr);
extern(C) int  swe_sol_eclipse_where(double tjd, int32 ifl, double *geopos, double *attr, char *serr);
extern(C) int  swe_lun_occult_where(double tjd, int32 ipl, const(char) *starname, int32 ifl, double *geopos, double *attr, char *serr);
extern(C) int  swe_sol_eclipse_how(double tjd, int32 ifl, double *geopos, double *attr, char *serr);
extern(C) int  swe_sol_eclipse_when_loc(double tjd_start, int32 ifl, double *geopos, double *tret, double *attr, int32 backward, char *serr);
extern(C) int  swe_lun_occult_when_loc(double tjd_start, int32 ipl, const(char) *starname, int32 ifl, double *geopos, double *tret, double *attr, int32 backward, char *serr);
extern(C) int  swe_sol_eclipse_when_glob(double tjd_start, int32 ifl, int32 ifltype, double *tret, int32 backward, char *serr);
extern(C) int  swe_lun_occult_when_glob(double tjd_start, int32 ipl, const(char) *starname, int32 ifl, int32 ifltype, double *tret, int32 backward, char *serr);
extern(C) int  swe_lun_eclipse_how(double tjd_ut, int32 ifl, double *geopos, double *attr, char *serr);
extern(C) int  swe_lun_eclipse_when(double tjd_start, int32 ifl, int32 ifltype, double *tret, int32 backward, char *serr);
extern(C) int  swe_lun_eclipse_when_loc(double tjd_start, int32 ifl, const(double) *geopos, double *tret, double *attr, int32 backward, char *serr);
extern(C) int  swe_pheno(double tjd, int32 ipl, int32 iflag, double *attr, char *serr);
extern(C) int  swe_pheno_ut(double tjd_ut, int32 ipl, int32 iflag, double *attr, char *serr);
extern(C) double swe_refrac(double inalt, double atpress, double attemp, int32 calc_flag);
extern(C) double swe_refrac_extended(double inalt, double geoalt, double atpress, double attemp, double lapse_rate, int32 calc_flag, double *dret);
extern(C) void swe_set_lapse_rate(double lapse_rate);
extern(C) void swe_azalt(double tjd_ut, int32 calc_flag, double *geopos, double atpress, double attemp, double *xin, double *xaz);
extern(C) void swe_azalt_rev(double tjd_ut, int32 calc_flag, double *geopos, double *xin, double *xout);
extern(C) int  swe_rise_trans_true_hor(double tjd_ut, int32 ipl, const(char) *starname, int32 epheflag, int32 rsmi,const(double) *geopos, double atpress, double attemp, double horhgt, double *tret, char *serr);
extern(C) int  swe_rise_trans(double tjd_ut, int32 ipl, const(char) *starname, int32 epheflag, int32 rsmi, const(double) *geopos, double atpress, double attemp, double *tret, char *serr);
extern(C) int  swe_nod_aps(double tjd_et, int32 ipl, int32 iflag, int32  method, double *xnasc, double *xndsc, double *xperi, double *xaphe, char *serr);
extern(C) int  swe_nod_aps_ut(double tjd_ut, int32 ipl, int32 iflag, int32  method, double *xnasc, double *xndsc, double *xperi, double *xaphe, char *serr);
extern(C) double  swe_deltat(double tjd);
extern(C) int  swe_time_equ(double tjd, double *te, char *serr);
extern(C) int  swe_lmt_to_lat(double tjd_lmt, double geolon, double *tjd_lat, char *serr);
extern(C) int  swe_lat_to_lmt(double tjd_lat, double geolon, double *tjd_lmt, char *serr);
extern(C) double  swe_sidtime0(double tjd_ut, double eps, double nut);
extern(C) double  swe_sidtime(double tjd_ut);
extern(C) void  swe_cotrans(double *xpo, double *xpn, double eps);
extern(C) void  swe_cotrans_sp(double *xpo, double *xpn, double eps);
extern(C) double  swe_get_tid_acc();
extern(C) void  swe_set_tid_acc(double t_acc);
extern(C) double  swe_degnorm(double x);
extern(C) double  swe_radnorm(double x);
extern(C) double  swe_rad_midp(double x1, double x0);
extern(C) double  swe_deg_midp(double x1, double x0);
extern(C) void  swe_split_deg(double ddeg, int32 roundflag, int32 *ideg, int32 *imin, int32 *isec, double *dsecfr, int32 *isgn);
extern(C) centisec  swe_csnorm(centisec p);
extern(C)  centisec  swe_difcsn (centisec p1, centisec p2);
extern(C) double  swe_difdegn (double p1, double p2);
extern(C) centisec swe_difcs2n(centisec p1, centisec p2);
extern(C) double  swe_difdeg2n(double p1, double p2);
extern(C) double  swe_difrad2n(double p1, double p2);
extern(C) centisec swe_csroundsec(centisec x);
extern(C) int  swe_d2l(double x);
extern(C) int  swe_day_of_week(double jd);
extern(C) char* swe_cs2timestr(CSEC t, int sep, AS_BOOL suppressZero, char *a);
extern(C) char* swe_cs2lonlatstr(CSEC t, char pchar, char mchar, char *s);
extern(C) char* swe_cs2degstr(CSEC t, char *a);

