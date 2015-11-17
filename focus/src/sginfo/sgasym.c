/*
  Space Group Info (c) 1994-97 Ralf W. Grosse-Kunstleve
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "sginfo.h"


static const char *AsymUnit1[] =
  {
    "0<=x<=1",
    "0<=y<=1",
    "0<=z<=1",
    NULL
  };
static const char *AsymUnit2[] =
  {
    "0<=x<=1/2",
    "0<=y<=1",
    "0<=z<=1",
    NULL
  };
static const char *AsymUnit3[] =
  {
    "0<=x<=1",
    "0<=y<=1",
    "0<=z<=1/2",
    NULL
  };
static const char *AsymUnit5[] =
  {
    "0<=x<=1/2",
    "0<=y<=1/2",
    "0<=z<=1",
    NULL
  };
static const char *AsymUnit6[] =
  {
    "0<=x<=1",
    "0<=y<=1/2",
    "0<=z<=1",
    NULL
  };
static const char *AsymUnit8[] =
  {
    "0<=x<=1",
    "0<=y<=1/4",
    "0<=z<=1",
    NULL
  };
static const char *AsymUnit12[] =
  {
    "0<=x<=1/2",
    "0<=y<=1/4",
    "0<=z<=1",
    NULL
  };
static const char *AsymUnit13[] =
  {
    "0<=x<=1/2",
    "0<=y<=1",
    "0<=z<=1/2",
    NULL
  };
static const char *AsymUnit15[] =
  {
    "0<=x<=1/2",
    "0<=y<=1/2",
    "0<=z<=1/2",
    NULL
  };
static const char *AsymUnit21[] =
  {
    "0<=x<=1/4",
    "0<=y<=1/2",
    "0<=z<=1",
    NULL
  };
static const char *AsymUnit22[] =
  {
    "0<=x<=1/4",
    "0<=y<=1/4",
    "0<=z<=1",
    NULL
  };
static const char *AsymUnit28[] =
  {
    "0<=x<=1/4",
    "0<=y<=1",
    "0<=z<=1",
    NULL
  };
static const char *AsymUnit46[] =
  {
    "0<=x<=1/4",
    "0<=y<=1",
    "0<=z<=1/2",
    NULL
  };
static const char *AsymUnit52[] =
  {
    "0<=x<=1",
    "0<=y<=1/4",
    "0<=z<=1/2",
    NULL
  };
static const char *AsymUnit53[] =
  {
    "0<=x<=1/2",
    "0<=y<=1",
    "0<=z<=1/4",
    NULL
  };
static const char *AsymUnit63[] =
  {
    "0<=x<=1/2",
    "0<=y<=1/2",
    "0<=z<=1/4",
    NULL
  };
static const char *AsymUnit64[] =
  {
    "0<=x<=1/4",
    "0<=y<=1/2",
    "0<=z<=1/2",
    NULL
  };
static const char *AsymUnit67[] =
  {
    "0<=x<=1/2",
    "0<=y<=1/4",
    "0<=z<=1/2",
    NULL
  };
static const char *AsymUnit69[] =
  {
    "0<=x<=1/4",
    "0<=y<=1/4",
    "0<=z<=1/2",
    NULL
  };
static const char *AsymUnit70[] =
  {
    "0<=x<=1/8",
    "0<=y<=1/4",
    "0<=z<=1",
    NULL
  };
static const char *AsymUnit91[] =
  {
    "0<=x<=1",
    "0<=y<=1",
    "0<=z<=1/8",
    NULL
  };
static const char *AsymUnit98[] =
  {
    "0<=x<=1/2",
    "0<=y<=1",
    "0<=z<=1/8",
    NULL
  };
static const char *AsymUnit99[] =
  {
    "0<=x<=1/2",
    "0<=y<=1/2",
    "0<=z<=1",
    "x<=y",
    NULL
  };
static const char *AsymUnit100[] =
  {
    "0<=x<=1/2",
    "0<=y<=1/2",
    "0<=z<=1",
    "y<=1/2-x",
    NULL
  };
static const char *AsymUnit107[] =
  {
    "0<=x<=1/2",
    "0<=y<=1/2",
    "0<=z<=1/2",
    "x<=y",
    NULL
  };
static const char *AsymUnit108[] =
  {
    "0<=x<=1/2",
    "0<=y<=1/2",
    "0<=z<=1/2",
    "y<=1/2-x",
    NULL
  };
static const char *AsymUnit134[] =
  {
    "0<=x<=1/2",
    "0<=y<=1",
    "0<=z<=1/4",
    "x<=y",
    "y<=1-x",
    NULL
  };
static const char *AsymUnit138[] =
  {
    "0<=x<=1/4",
    "0<=y<=1/2",
    "0<=z<=1",
    "x<=y",
    "y<=1/2-x",
    NULL
  };
static const char *AsymUnit139[] =
  {
    "0<=x<=1/2",
    "0<=y<=1/2",
    "0<=z<=1/4",
    "x<=y",
    NULL
  };
static const char *AsymUnit140[] =
  {
    "0<=x<=1/2",
    "0<=y<=1/2",
    "0<=z<=1/4",
    "y<=1/2-x",
    NULL
  };
static const char *AsymUnit141[] =
  {
    "0<=x<=1/2",
    "0<=y<=1/2",
    "0<=z<=1/8",
    NULL
  };
static const char *AsymUnit143[] =
  {
    "0<=x<=2/3",
    "0<=y<=2/3",
    "0<=z<=1",
    "x<=(1+y)/2",
    "y<=min(1-x,(1+x)/2)",
    NULL
  };
static const char *AsymUnit144[] =
  {
    "0<=x<=1",
    "0<=y<=1",
    "0<=z<=1/3",
    NULL
  };
static const char *AsymUnit146[] =
  {
    "0<=x<=2/3",
    "0<=y<=2/3",
    "0<=z<=1/3",
    "x<=(1+y)/2",
    "y<=min(1-x,(1+x)/2)",
    NULL
  };
static const char *AsymUnit147[] =
  {
    "0<=x<=2/3",
    "0<=y<=2/3",
    "0<=z<=1/2",
    "x<=(1+y)/2",
    "y<=min(1-x,(1+x)/2)",
    NULL
  };
static const char *AsymUnit148[] =
  {
    "0<=x<=2/3",
    "0<=y<=2/3",
    "0<=z<=1/6",
    "x<=(1+y)/2",
    "y<=min(1-x,(1+x)/2)",
    NULL
  };
static const char *AsymUnit151[] =
  {
    "0<=x<=1",
    "0<=y<=1",
    "0<=z<=1/6",
    NULL
  };
static const char *AsymUnit156[] =
  {
    "0<=x<=2/3",
    "0<=y<=2/3",
    "0<=z<=1",
    "x<=2 y",
    "y<=min(1-x,2 x)",
    NULL
  };
static const char *AsymUnit157[] =
  {
    "0<=x<=2/3",
    "0<=y<=1/2",
    "0<=z<=1",
    "x<=(y+1)/2",
    "y<=min(1-x,x)",
    NULL
  };
static const char *AsymUnit160[] =
  {
    "0<=x<=2/3",
    "0<=y<=2/3",
    "0<=z<=1/3",
    "x<=2 y",
    "y<=min(1-x,2 x)",
    NULL
  };
static const char *AsymUnit162[] =
  {
    "0<=x<=2/3",
    "0<=y<=1/2",
    "0<=z<=1/2",
    "x<=(1+y)/2",
    "y<=min(1-x,x)",
    NULL
  };
static const char *AsymUnit163[] =
  {
    "0<=x<=2/3",
    "0<=y<=2/3",
    "0<=z<=1/4",
    "x<=(1+y)/2",
    "y<=min(1-x,(1+x)/2)",
    NULL
  };
static const char *AsymUnit164[] =
  {
    "0<=x<=2/3",
    "0<=y<=1/3",
    "0<=z<=1",
    "x<=(1+y)/2",
    "y<=x/2",
    NULL
  };
static const char *AsymUnit166[] =
  {
    "0<=x<=2/3",
    "0<=y<=2/3",
    "0<=z<=1/6",
    "x<=2 y",
    "y<=min(1-x,2 x)",
    NULL
  };
static const char *AsymUnit167[] =
  {
    "0<=x<=2/3",
    "0<=y<=2/3",
    "0<=z<=1/12",
    "x<=(1+y)/2",
    "y<=min(1-x,(1+x)/2)",
    NULL
  };
static const char *AsymUnit168[] =
  {
    "0<=x<=2/3",
    "0<=y<=1/2",
    "0<=z<=1",
    "x<=(1+y)/2",
    "y<=min(1-x,x)",
    NULL
  };
static const char *AsymUnit171[] =
  {
    "0<=x<=1",
    "0<=y<=1",
    "0<=z<=1/3",
    "y<=x",
    NULL
  };
static const char *AsymUnit178[] =
  {
    "0<=x<=1",
    "0<=y<=1",
    "0<=z<=1/12",
    NULL
  };
static const char *AsymUnit180[] =
  {
    "0<=x<=1",
    "0<=y<=1",
    "0<=z<=1/6",
    "y<=x",
    NULL
  };
static const char *AsymUnit187[] =
  {
    "0<=x<=2/3",
    "0<=y<=2/3",
    "0<=z<=1/2",
    "x<=2 y",
    "y<=min(1-x,2 x)",
    NULL
  };
static const char *AsymUnit191[] =
  {
    "0<=x<=2/3",
    "0<=y<=1/3",
    "0<=z<=1/2",
    "x<=(1+y)/2",
    "y<=x/2",
    NULL
  };
static const char *AsymUnit192[] =
  {
    "0<=x<=2/3",
    "0<=y<=1/2",
    "0<=z<=1/4",
    "x<=(1+y)/2",
    "y<=min(1-x,x)",
    NULL
  };
static const char *AsymUnit194[] =
  {
    "0<=x<=2/3",
    "0<=y<=2/3",
    "0<=z<=1/4",
    "x<=2 y",
    "y<=min(1-x,2 x)",
    NULL
  };
static const char *AsymUnit195[] =
  {
    "0<=x<=1",
    "0<=y<=1",
    "0<=z<=1/2",
    "y<=1-x",
    "z<=min(x,y)",
    NULL
  };
static const char *AsymUnit196[] =
  {
    "0<=x<=1/2",
    "0<=y<=1/2",
    "-1/4<=z<=1/4",
    "y<=x",
    "max(x-1/2,-y)<=z<=min(1/2-x,y)",
    NULL
  };
static const char *AsymUnit197[] =
  {
    "0<=x<=1",
    "0<=y<=1/2",
    "0<=z<=1/2",
    "y<=min(x,1-x)",
    "z<=y",
    NULL
  };
static const char *AsymUnit198[] =
  {
    "0<=x<=1/2",
    "0<=y<=1/2",
    "-1/2<=z<=1/2",
    "max(x-1/2,-y)<=z<=min(x,y)",
    NULL
  };
static const char *AsymUnit199[] =
  {
    "0<=x<=1/2",
    "0<=y<=1/2",
    "0<=z<=1/2",
    "z<=min(x,y)",
    NULL
  };
static const char *AsymUnit202[] =
  {
    "0<=x<=1/2",
    "0<=y<=1/2",
    "0<=z<=1/4",
    "y<=x",
    "z<=min(1/2-x,y)",
    NULL
  };
static const char *AsymUnit203[] =
  {
    "0<=x<=1/2",
    "0<=y<=1/4",
    "-1/4<=z<=1/4",
    "y<=min(x,1/2-x)",
    "-y<=z<=y",
    NULL
  };
static const char *AsymUnit204[] =
  {
    "0<=x<=1/2",
    "0<=y<=1/2",
    "0<=z<=1/2",
    "y<=x",
    "z<=y",
    NULL
  };
static const char *AsymUnit206[] =
  {
    "0<=x<=1/2",
    "0<=y<=1/2",
    "0<=z<=1/4",
    "z<=min(x,1/2-x,1/2-y)",
    NULL
  };
static const char *AsymUnit208[] =
  {
    "0<=x<=1/2",
    "0<=y<=1/2",
    "-1/4<=z<=1/4",
    "max(-x,x-1/2,-y,y-1/2)<=z<=min(x,1/2-x,y,1/2-y)",
    NULL
  };
static const char *AsymUnit210[] =
  {
    "0<=x<=1/2",
    "-1/8<=y<=1/8",
    "-1/8<=z<=1/8",
    "y<=min(x,1/2-x)",
    "-y<=z<=min(x,1/2-x)",
    NULL
  };
static const char *AsymUnit211[] =
  {
    "0<=x<=1/2",
    "0<=y<=1/2",
    "0<=z<=1/4",
    "z<=min(x,1/2-x,y,1/2-y)",
    NULL
  };
static const char *AsymUnit212[] =
  {
    "0<=x<=1/2",
    "0<=y<=3/4",
    "-1/2<=z<=1/4",
    "max(-y,x-1/2)<=z<=min(1/2-y,2 x-y,2 y-x,y-2 x+1/2)",
    NULL
  };
static const char *AsymUnit213[] =
  {
    "-1/4<=x<=1/2",
    "0<=y<=3/4",
    "0<=z<=1/2",
    "x<=y<=x+1/2",
    "(y-x)/2<=z<=min(y,(-4 x-2 y+3)/2,(3-2 x-2 y)/4)",
    NULL
  };
static const char *AsymUnit214[] =
  {
    "-3/8<=x<=1/8",
    "-1/8<=y<=1/8",
    "-1/8<=z<=3/8",
    "max(x,y,y-x-1/8)<=z<=y+1/4",
    NULL
  };
static const char *AsymUnit220[] =
  {
    "1/4<=x<=1/2",
    "1/4<=y<=1/2",
    "0<=z<=1/2",
    "z<=min(x,y)",
    NULL
  };
static const char *AsymUnit225[] =
  {
    "0<=x<=1/2",
    "0<=y<=1/4",
    "0<=z<=1/4",
    "y<=min(x,1/2-x)",
    "z<=y",
    NULL
  };
static const char *AsymUnit227[] =
  {
    "0<=x<=1/2",
    "0<=y<=1/8",
    "-1/8<=z<=1/8",
    "y<=min(x,1/2-x)",
    "-y<=z<=y",
    NULL
  };
static const char *AsymUnit230[] =
  {
    "-1/8<=x<=1/8",
    "-1/8<=y<=1/8",
    "0<=z<=1/4",
    "max(x,-x,y,-y)<=z",
    NULL
  };


static const char **AsymUnits[] =
  {
    NULL,
    /*   1 */ AsymUnit1,
    /*   2 */ AsymUnit2,
    /*   3 */ AsymUnit3,
    /*   4 */ AsymUnit3,
    /*   5 */ AsymUnit5,
    /*   6 */ AsymUnit6,
    /*   7 */ AsymUnit6,
    /*   8 */ AsymUnit8,
    /*   9 */ AsymUnit8,
    /*  10 */ AsymUnit5,
    /*  11 */ AsymUnit8,
    /*  12 */ AsymUnit12,
    /*  13 */ AsymUnit13,
    /*  14 */ AsymUnit8,
    /*  15 */ AsymUnit15,
    /*  16 */ AsymUnit5,
    /*  17 */ AsymUnit5,
    /*  18 */ AsymUnit5,
    /*  19 */ AsymUnit5,
    /*  20 */ AsymUnit15,
    /*  21 */ AsymUnit21,
    /*  22 */ AsymUnit22,
    /*  23 */ AsymUnit15,
    /*  24 */ AsymUnit15,
    /*  25 */ AsymUnit5,
    /*  26 */ AsymUnit5,
    /*  27 */ AsymUnit5,
    /*  28 */ AsymUnit28,
    /*  29 */ AsymUnit28,
    /*  30 */ AsymUnit13,
    /*  31 */ AsymUnit5,
    /*  32 */ AsymUnit5,
    /*  33 */ AsymUnit5,
    /*  34 */ AsymUnit5,
    /*  35 */ AsymUnit21,
    /*  36 */ AsymUnit15,
    /*  37 */ AsymUnit21,
    /*  38 */ AsymUnit15,
    /*  39 */ AsymUnit12,
    /*  40 */ AsymUnit21,
    /*  41 */ AsymUnit15,
    /*  42 */ AsymUnit22,
    /*  43 */ AsymUnit22,
    /*  44 */ AsymUnit15,
    /*  45 */ AsymUnit15,
    /*  46 */ AsymUnit46,
    /*  47 */ AsymUnit15,
    /*  48 */ AsymUnit21,
    /*  49 */ AsymUnit15,
    /*  50 */ AsymUnit15,
    /*  51 */ AsymUnit21,
    /*  52 */ AsymUnit52,
    /*  53 */ AsymUnit53,
    /*  54 */ AsymUnit15,
    /*  55 */ AsymUnit15,
    /*  56 */ AsymUnit46,
    /*  57 */ AsymUnit53,
    /*  58 */ AsymUnit15,
    /*  59 */ AsymUnit15,
    /*  60 */ AsymUnit15,
    /*  61 */ AsymUnit15,
    /*  62 */ AsymUnit12,
    /*  63 */ AsymUnit63,
    /*  64 */ AsymUnit64,
    /*  65 */ AsymUnit64,
    /*  66 */ AsymUnit64,
    /*  67 */ AsymUnit67,
    /*  68 */ AsymUnit64,
    /*  69 */ AsymUnit69,
    /*  70 */ AsymUnit70,
    /*  71 */ AsymUnit64,
    /*  72 */ AsymUnit64,
    /*  73 */ AsymUnit64,
    /*  74 */ AsymUnit22,
    /*  75 */ AsymUnit5,
    /*  76 */ AsymUnit5,
    /*  77 */ AsymUnit5,
    /*  78 */ AsymUnit5,
    /*  79 */ AsymUnit15,
    /*  80 */ AsymUnit53,
    /*  81 */ AsymUnit5,
    /*  82 */ AsymUnit15,
    /*  83 */ AsymUnit15,
    /*  84 */ AsymUnit15,
    /*  85 */ AsymUnit15,
    /*  86 */ AsymUnit53,
    /*  87 */ AsymUnit63,
    /*  88 */ AsymUnit22,
    /*  89 */ AsymUnit15,
    /*  90 */ AsymUnit15,
    /*  91 */ AsymUnit91,
    /*  92 */ AsymUnit91,
    /*  93 */ AsymUnit53,
    /*  94 */ AsymUnit15,
    /*  95 */ AsymUnit91,
    /*  96 */ AsymUnit91,
    /*  97 */ AsymUnit63,
    /*  98 */ AsymUnit98,
    /*  99 */ AsymUnit99,
    /* 100 */ AsymUnit100,
    /* 101 */ AsymUnit99,
    /* 102 */ AsymUnit99,
    /* 103 */ AsymUnit15,
    /* 104 */ AsymUnit15,
    /* 105 */ AsymUnit15,
    /* 106 */ AsymUnit15,
    /* 107 */ AsymUnit107,
    /* 108 */ AsymUnit108,
    /* 109 */ AsymUnit63,
    /* 110 */ AsymUnit63,
    /* 111 */ AsymUnit99,
    /* 112 */ AsymUnit15,
    /* 113 */ AsymUnit100,
    /* 114 */ AsymUnit15,
    /* 115 */ AsymUnit15,
    /* 116 */ AsymUnit53,
    /* 117 */ AsymUnit15,
    /* 118 */ AsymUnit53,
    /* 119 */ AsymUnit63,
    /* 120 */ AsymUnit63,
    /* 121 */ AsymUnit107,
    /* 122 */ AsymUnit98,
    /* 123 */ AsymUnit107,
    /* 124 */ AsymUnit63,
    /* 125 */ AsymUnit108,
    /* 126 */ AsymUnit63,
    /* 127 */ AsymUnit108,
    /* 128 */ AsymUnit63,
    /* 129 */ AsymUnit108,
    /* 130 */ AsymUnit63,
    /* 131 */ AsymUnit63,
    /* 132 */ AsymUnit107,
    /* 133 */ AsymUnit63,
    /* 134 */ AsymUnit134,
    /* 135 */ AsymUnit63,
    /* 136 */ AsymUnit107,
    /* 137 */ AsymUnit63,
    /* 138 */ AsymUnit138,
    /* 139 */ AsymUnit139,
    /* 140 */ AsymUnit140,
    /* 141 */ AsymUnit141,
    /* 142 */ AsymUnit141,
    /* 143 */ AsymUnit143,
    /* 144 */ AsymUnit144,
    /* 145 */ AsymUnit144,
    /* 146 */ AsymUnit146,
    /* 147 */ AsymUnit147,
    /* 148 */ AsymUnit148,
    /* 149 */ AsymUnit147,
    /* 150 */ AsymUnit147,
    /* 151 */ AsymUnit151,
    /* 152 */ AsymUnit151,
    /* 153 */ AsymUnit151,
    /* 154 */ AsymUnit151,
    /* 155 */ AsymUnit148,
    /* 156 */ AsymUnit156,
    /* 157 */ AsymUnit157,
    /* 158 */ AsymUnit147,
    /* 159 */ AsymUnit147,
    /* 160 */ AsymUnit160,
    /* 161 */ AsymUnit148,
    /* 162 */ AsymUnit162,
    /* 163 */ AsymUnit163,
    /* 164 */ AsymUnit164,
    /* 165 */ AsymUnit163,
    /* 166 */ AsymUnit166,
    /* 167 */ AsymUnit167,
    /* 168 */ AsymUnit168,
    /* 169 */ AsymUnit151,
    /* 170 */ AsymUnit151,
    /* 171 */ AsymUnit171,
    /* 172 */ AsymUnit171,
    /* 173 */ AsymUnit147,
    /* 174 */ AsymUnit147,
    /* 175 */ AsymUnit162,
    /* 176 */ AsymUnit163,
    /* 177 */ AsymUnit162,
    /* 178 */ AsymUnit178,
    /* 179 */ AsymUnit178,
    /* 180 */ AsymUnit180,
    /* 181 */ AsymUnit180,
    /* 182 */ AsymUnit163,
    /* 183 */ AsymUnit164,
    /* 184 */ AsymUnit162,
    /* 185 */ AsymUnit162,
    /* 186 */ AsymUnit164,
    /* 187 */ AsymUnit187,
    /* 188 */ AsymUnit163,
    /* 189 */ AsymUnit162,
    /* 190 */ AsymUnit163,
    /* 191 */ AsymUnit191,
    /* 192 */ AsymUnit192,
    /* 193 */ AsymUnit192,
    /* 194 */ AsymUnit194,
    /* 195 */ AsymUnit195,
    /* 196 */ AsymUnit196,
    /* 197 */ AsymUnit197,
    /* 198 */ AsymUnit198,
    /* 199 */ AsymUnit199,
    /* 200 */ AsymUnit199,
    /* 201 */ AsymUnit197,
    /* 202 */ AsymUnit202,
    /* 203 */ AsymUnit203,
    /* 204 */ AsymUnit204,
    /* 205 */ AsymUnit199,
    /* 206 */ AsymUnit206,
    /* 207 */ AsymUnit197,
    /* 208 */ AsymUnit208,
    /* 209 */ AsymUnit203,
    /* 210 */ AsymUnit210,
    /* 211 */ AsymUnit211,
    /* 212 */ AsymUnit212,
    /* 213 */ AsymUnit213,
    /* 214 */ AsymUnit214,
    /* 215 */ AsymUnit197,
    /* 216 */ AsymUnit203,
    /* 217 */ AsymUnit204,
    /* 218 */ AsymUnit199,
    /* 219 */ AsymUnit203,
    /* 220 */ AsymUnit220,
    /* 221 */ AsymUnit204,
    /* 222 */ AsymUnit204,
    /* 223 */ AsymUnit211,
    /* 224 */ AsymUnit196,
    /* 225 */ AsymUnit225,
    /* 226 */ AsymUnit225,
    /* 227 */ AsymUnit227,
    /* 228 */ AsymUnit227,
    /* 229 */ AsymUnit202,
    /* 230 */ AsymUnit230
  };


const char **StdAsymUnit(int SgNumber, const char *HallSymbol)
{
  const T_TabSgName  *TSgN;


  if (SgNumber < 1 || SgNumber > 230) return NULL;

  if (HallSymbol == NULL)
    return AsymUnits[SgNumber];

      TSgN = FindTabSgNameEntry(SgNumber, NULL, 'A');
  if (TSgN == NULL) {
    SetSgError("Internal Error: Corrupt FindTabSgNameEntry()");
    return NULL;
  }

  if (strcmp(HallSymbol, TSgN->HallSymbol) == 0)
    return AsymUnits[SgNumber];

  return NULL;
}
