#include "energy_const.h"
#include <string>

using namespace std;

static int hairpin37[31] = {
  INF, INF, INF, 570, 560, 560, 540, 590, 560, 640, 650,
       660, 670, 678, 686, 694, 701, 707, 713, 719, 725,
       730, 735, 740, 744, 749, 753, 757, 761, 765, 769};

static string Tetraloops = 
  "GGGGAC GGUGAC CGAAAG GGAGAC CGCAAG GGAAAC CGGAAG CUUCGG CGUGAG CGAAGG CUACGG GGCAAC CGCGAG UGAGAG CGAGAG AGAAAU CGUAAG CUAACG UGAAAG GGAAGC GGGAAC UGAAAA AGCAAU AGUAAU CGGGAG AGUGAU GGCGAC GGGAGC GUGAAC UGGAAA ";

static int tetra_energy37[200] = {
  -300, -300, -300, -300, -300, -300,
  -300, -300, -300, -250, -250, -250,
  -250, -250, -200, -200, -200, -200,
  -200, -150, -150, -150, -150, -150,
  -150, -150, -150, -150, -150, -150};

static int mismatchH37[7][5][5] =
{ /* @@ */
  {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
  { /* CG */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   { -90, -150, -150, -140, -180}, /* A@  AA  AC  AG  AU */
   { -90, -100,  -90, -290,  -80}, /* C@  CA  CC  CG  CU */
   { -90, -220, -200, -160, -110}, /* G@  GA  GC  GG  GU */
   { -90, -170, -140, -180, -200}},/* U@  UA  UC  UG  UU */
  { /* GC */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   { -70, -110, -150, -130, -210}, /* A@  AA  AC  AG  AU */
   { -70, -110,  -70, -240,  -50}, /* C@  CA  CC  CG  CU */
   { -70, -240, -290, -140, -120}, /* G@  GA  GC  GG  GU */
   { -70, -190, -100, -220, -150}},/* U@  UA  UC  UG  UU */
  { /* GU */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,   20,  -50,  -30,  -30}, /* A@  AA  AC  AG  AU */
   {   0,  -10,  -20, -150,  -20}, /* C@  CA  CC  CG  CU */
   {   0,  -90, -110,  -30,    0}, /* G@  GA  GC  GG  GU */
   {   0,  -30,  -30,  -40, -110}},/* U@  UA  UC  UG  UU */
  { /* UG */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,  -50,  -30,  -60,  -50}, /* A@  AA  AC  AG  AU */
   {   0,  -20,  -10, -170,    0}, /* C@  CA  CC  CG  CU */
   {   0,  -80, -120,  -30,  -70}, /* G@  GA  GC  GG  GU */
   {   0,  -60,  -10,  -60,  -80}},/* U@  UA  UC  UG  UU */
  { /* AU */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,  -30,  -50,  -30,  -30}, /* A@  AA  AC  AG  AU */
   {   0,  -10,  -20, -150,  -20}, /* C@  CA  CC  CG  CU */
   {   0, -110, -120,  -20,   20}, /* G@  GA  GC  GG  GU */
   {   0,  -30,  -30,  -60, -110}},/* U@  UA  UC  UG  UU */
  { /* UA */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,  -50,  -30,  -60,  -50}, /* A@  AA  AC  AG  AU */
   {   0,  -20,  -10, -120,   -0}, /* C@  CA  CC  CG  CU */
   {   0, -140, -120,  -70,  -20}, /* G@  GA  GC  GG  GU */
   {   0,  -30,  -10,  -50,  -80}}/* U@  UA  UC  UG  UU */
};

static int stack37[7][7] =
/*          CG     GC     GU     UG     AU     UA  */
{ {  INF,   INF,   INF,   INF,   INF,   INF,   INF},
  {  INF,  -240,  -330,  -210,  -140,  -210,  -210},
  {  INF,  -330,  -340,  -250,  -150,  -220,  -240},
  {  INF,  -210,  -250,   130,   -50,  -140,  -130},
  {  INF,  -140,  -150,   -50,    30,   -60,  -100},
  {  INF,  -210,  -220,  -140,   -60,  -110,   -90},
  {  INF,  -210,  -240,  -130,  -100,   -90,  -130} };

static int bulge37[31] = {
  INF, 380, 280, 320, 360, 400, 440, 459, 470, 480, 490,
       500, 510, 519, 527, 534, 541, 548, 554, 560, 565,
	   571, 576, 580, 585, 589, 594, 598, 602, 605, 609};

static int TerminalAU = 50;

static int internal_loop37[31] = {
  INF, INF, 410, 510, 170, 180, 200, 220, 230, 240, 250,
       260, 270, 278, 286, 294, 301, 307, 313, 319, 325,
       330, 335, 340, 345, 349, 353, 357, 361, 365, 369};

static int mismatchI37[7][5][5] =
{ /* @@ */
  {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
  { /* CG */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,    0,    0, -110,    0}, /* A@  AA  AC  AG  AU */
   {   0,    0,    0,    0,    0}, /* C@  CA  CC  CG  CU */
   {   0, -110,    0,    0,    0}, /* G@  GA  GC  GG  GU */
   {   0,    0,    0,    0,  -70}},/* U@  UA  UC  UG  UU */
  { /* GC */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,    0,    0, -110,    0}, /* A@  AA  AC  AG  AU */
   {   0,    0,    0,    0,    0}, /* C@  CA  CC  CG  CU */
   {   0, -110,    0,    0,    0}, /* G@  GA  GC  GG  GU */
   {   0,    0,    0,    0,  -70}},/* U@  UA  UC  UG  UU */
  { /* GU */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,   70,   70,  -40,   70}, /* A@  AA  AC  AG  AU */
   {   0,   70,   70,   70,   70}, /* C@  CA  CC  CG  CU */
   {   0,  -40,   70,   70,   70}, /* G@  GA  GC  GG  GU */
   {   0,   70,   70,   70,    0}},/* U@  UA  UC  UG  UU */
  { /* UG */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,   70,   70,  -40,   70}, /* A@  AA  AC  AG  AU */
   {   0,   70,   70,   70,   70}, /* C@  CA  CC  CG  CU */
   {   0,  -40,   70,   70,   70}, /* G@  GA  GC  GG  GU */
   {   0,   70,   70,   70,    0}},/* U@  UA  UC  UG  UU */
  { /* AU */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,   70,   70,  -40,   70}, /* A@  AA  AC  AG  AU */
   {   0,   70,   70,   70,   70}, /* C@  CA  CC  CG  CU */
   {   0,  -40,   70,   70,   70}, /* G@  GA  GC  GG  GU */
   {   0,   70,   70,   70,    0}},/* U@  UA  UC  UG  UU */
  { /* UA */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,   70,   70,  -40,   70}, /* A@  AA  AC  AG  AU */
   {   0,   70,   70,   70,   70}, /* C@  CA  CC  CG  CU */
   {   0,  -40,   70,   70,   70}, /* G@  GA  GC  GG  GU */
   {   0,   70,   70,   70,    0}},/* U@  UA  UC  UG  UU */
};

static int ML_closing37 = 340;
static int ML_intern37 = 40;

static int dangle5_37[8][5]=
{/*   @     A     C     G     U   */
   { INF,  INF,  INF,  INF,  INF}, /* no pair */
   { INF,  -50,  -30,  -20,  -10}, /* CG  (stacks on C) */
   { INF,  -20,  -30,   -0,   -0}, /* GC  (stacks on G) */
   { INF,  -30,  -30,  -40,  -20}, /* GU */
   { INF,  -30,  -10,  -20,  -20}, /* UG */
   { INF,  -30,  -30,  -40,  -20}, /* AU */
   { INF,  -30,  -10,  -20,  -20}, /* UA */
   {   0,    0,     0,    0,   0}  /*  @ */
};

/* 3' dangling ends (unpaired base stacks on second paired base */
static int dangle3_37[8][5]=
{/*   @     A     C     G     U   */
   { INF,  INF,  INF,  INF,  INF},  /* no pair */
   { INF, -110,  -40, -130,  -60},  /* CG  (stacks on G) */
   { INF, -170,  -80, -170, -120},  /* GC */
   { INF,  -70,  -10,  -70,  -10},  /* GU */
   { INF,  -80,  -50,  -80,  -60},  /* UG */
   { INF,  -70,  -10,  -70,  -10},  /* AU */
   { INF,  -80,  -50,  -80,  -60},  /* UA */
   {   0,    0,     0,    0,   0}   /*  @ */
};

static int ML_BASE37 = 0;

static int MAX_NINIO = 300;
static int F_ninio37 = 50;

#include "intloops.h"