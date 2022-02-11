(* Content-type: application/vnd.wolfram.cdf.text *)

(*** Wolfram CDF File ***)
(* http://www.wolfram.com/cdf *)

(* CreatedBy='Mathematica 9.0' *)

(*************************************************************************)
(*                                                                       *)
(*  The Mathematica License under which this file was created prohibits  *)
(*  restricting third parties in receipt of this file from republishing  *)
(*  or redistributing it by any means, including but not limited to      *)
(*  rights management or terms of use, without the express consent of    *)
(*  Wolfram Research, Inc. For additional information concerning CDF     *)
(*  licensing and redistribution see:                                    *)
(*                                                                       *)
(*        www.wolfram.com/cdf/adopting-cdf/licensing-options.html        *)
(*                                                                       *)
(*************************************************************************)

(*CacheID: 234*)
(* Internal cache information:
NotebookFileLineBreakTest
NotebookFileLineBreakTest
NotebookDataPosition[      1063,         20]
NotebookDataLength[     83085,       1810]
NotebookOptionsPosition[     82402,       1765]
NotebookOutlinePosition[     82824,       1781]
CellTagsIndexPosition[     82781,       1778]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell["Define functions for MAD data", "Title",
 Background->RGBColor[1, 1, 0]],

Cell[CellGroupData[{

Cell["\<\
Extract the coordinates: X,Y from survey data base\
\>", "Subsection",
 CellChangeTimes->{{3.6053446581155024`*^9, 3.6053446618751745`*^9}, {
  3.6056377946897187`*^9, 3.6056377948457217`*^9}, {3.690263978941309*^9, 
  3.690263987053465*^9}},
 Background->RGBColor[0.803922, 0.682353, 0.196078]],

Cell[BoxData[{
 RowBox[{
  RowBox[{"Clear", "[", "GetXYDB", "]"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"GetXYDB", "[", "x_", "]"}], ":=", 
  RowBox[{"Module", "[", "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"{", "f", "}"}], ",", "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{
      RowBox[{"f", "[", "y_", "]"}], ":=", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"y", "[", 
         RowBox[{"[", "2", "]"}], "]"}], ",", 
        RowBox[{"y", "[", 
         RowBox[{"[", "3", "]"}], "]"}]}], "}"}]}], ";", 
     "\[IndentingNewLine]", 
     RowBox[{"f", "/@", "x"}]}]}], "\[IndentingNewLine]", "]"}]}]}], "Input",
 CellChangeTimes->{{3.6053446798623204`*^9, 3.6053447440575547`*^9}, {
  3.6053467519438696`*^9, 3.605346770633029*^9}, {3.605637117418295*^9, 
  3.605637127090481*^9}, {3.605637249989644*^9, 3.6056372501612473`*^9}, {
  3.605637294434899*^9, 3.6056373173205385`*^9}, {3.6056378011794434`*^9, 
  3.605637811694046*^9}, {3.6902639342464495`*^9, 3.6902639531228123`*^9}}]
}, Open  ]],

Cell[CellGroupData[{

Cell["\<\
Extract the coordinates: X,Y from MADX\
\>", "Subsection",
 CellChangeTimes->{{3.6053446581155024`*^9, 3.6053446618751745`*^9}, {
  3.6056377946897187`*^9, 3.6056377948457217`*^9}, {3.6862916403060255`*^9, 
  3.6862916419460583`*^9}, {3.6902639899395204`*^9, 3.6902639957584324`*^9}},
 Background->RGBColor[0.803922, 0.682353, 0.196078]],

Cell[BoxData[{
 RowBox[{
  RowBox[{"Clear", "[", "GetXYMADX", "]"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"GetXYMADX", "[", "x_", "]"}], ":=", 
  RowBox[{"Module", "[", "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"{", "f", "}"}], ",", "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{
      RowBox[{"f", "[", "y_", "]"}], ":=", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"y", "[", 
         RowBox[{"[", "4", "]"}], "]"}], ",", 
        RowBox[{"y", "[", 
         RowBox[{"[", "5", "]"}], "]"}]}], "}"}]}], ";", 
     "\[IndentingNewLine]", 
     RowBox[{"f", "/@", "x"}]}]}], "\[IndentingNewLine]", "]"}]}]}], "Input",
 CellChangeTimes->{{3.6053446798623204`*^9, 3.6053447440575547`*^9}, {
  3.6053467519438696`*^9, 3.605346770633029*^9}, {3.605637117418295*^9, 
  3.605637127090481*^9}, {3.605637249989644*^9, 3.6056372501612473`*^9}, {
  3.605637294434899*^9, 3.6056373173205385`*^9}, {3.6056378011794434`*^9, 
  3.605637811694046*^9}, {3.686291645456128*^9, 3.6862916547163134`*^9}, {
  3.690264001733347*^9, 3.690264013683177*^9}, {3.6902644054379106`*^9, 
  3.6902644080899615`*^9}}]
}, Open  ]],

Cell[CellGroupData[{

Cell["PlotExplorer", "Subsection",
 CellChangeTimes->{{3.6053446581155024`*^9, 3.6053446618751745`*^9}, {
  3.6056377946897187`*^9, 3.6056377948457217`*^9}, {3.6862916403060255`*^9, 
  3.6862916419460583`*^9}, {3.6902639899395204`*^9, 3.6902639957584324`*^9}, {
  3.6902682058731995`*^9, 3.690268209193266*^9}},
 Background->RGBColor[0.803922, 0.682353, 0.196078]],

Cell[BoxData[{
 RowBox[{
  RowBox[{
   RowBox[{"Attributes", "[", "PlotExplorer", "]"}], "=", 
   RowBox[{"{", "HoldFirst", "}"}]}], ";"}], "\n", 
 RowBox[{
  RowBox[{
   RowBox[{"PlotExplorer", "[", "defPlot_", "]"}], ":=", 
   RowBox[{"DynamicModule", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{
      "held", ",", "head", ",", "plot", ",", "range", ",", "zoomRange", ",", 
       "defRange", ",", "defSize", ",", "defEpilog", ",", "size", ",", "pt1", 
       ",", "pt2", ",", "spt1", ",", "dist", ",", 
       RowBox[{"rect", "=", 
        RowBox[{"{", "}"}]}], ",", 
       RowBox[{"dragged", "=", "False"}], ",", 
       RowBox[{"overButton", "=", "False"}], ",", 
       RowBox[{"pushedButton", "=", "False"}], ",", 
       RowBox[{"overHor", "=", "False"}], ",", 
       RowBox[{"overVer", "=", "False"}], ",", 
       RowBox[{"activeHor", "=", "False"}], ",", 
       RowBox[{"activeVer", "=", "False"}], ",", 
       RowBox[{"xPos", "=", ".5"}], ",", 
       RowBox[{"yPos", "=", ".5"}], ",", 
       RowBox[{"minDistance", "=", "0.1"}], ",", 
       RowBox[{"sliderRatio", "=", "0.05"}], ",", 
       RowBox[{"panFactor", "=", "10"}], ",", "reset", ",", "distance", ",", 
       "slider", ",", "replot", ",", 
       RowBox[{"coord", "=", "None"}]}], "}"}], ",", 
     RowBox[{
      RowBox[{"held", "=", 
       RowBox[{"Hold", "@", "defPlot"}]}], ";", "\[IndentingNewLine]", 
      RowBox[{"head", "=", 
       RowBox[{"Part", "[", 
        RowBox[{"held", ",", "1", ",", "0"}], "]"}]}], ";", 
      "\[IndentingNewLine]", 
      RowBox[{"plot", "=", "defPlot"}], ";", "\[IndentingNewLine]", 
      RowBox[{
       RowBox[{"{", 
        RowBox[{"defEpilog", ",", "defRange"}], "}"}], "=", 
       RowBox[{
        RowBox[{"{", 
         RowBox[{"Epilog", ",", "PlotRange"}], "}"}], "/.", 
        RowBox[{"Quiet", "@", 
         RowBox[{"AbsoluteOptions", "@", "plot"}]}]}]}], ";", 
      "\[IndentingNewLine]", 
      RowBox[{"defEpilog", "=", 
       RowBox[{"defEpilog", "/.", 
        RowBox[{"Epilog", "\[Rule]", 
         RowBox[{"{", "}"}]}]}]}], ";", "\[IndentingNewLine]", 
      RowBox[{"size", "=", 
       RowBox[{"defSize", "=", 
        RowBox[{"Rasterize", "[", 
         RowBox[{"plot", ",", "\"\<RasterSize\>\""}], "]"}]}]}], ";", 
      "\[IndentingNewLine]", 
      RowBox[{"range", "=", 
       RowBox[{"zoomRange", "=", "defRange"}]}], ";", "\[IndentingNewLine]", 
      RowBox[{"Deploy", "@", 
       RowBox[{"EventHandler", "[", 
        RowBox[{
         RowBox[{"Dynamic", "[", 
          RowBox[{
           RowBox[{"MouseAppearance", "[", 
            RowBox[{
             RowBox[{"Show", "[", 
              RowBox[{"plot", ",", 
               RowBox[{"PlotRange", "\[Rule]", 
                RowBox[{"Dynamic", "@", "range"}]}], ",", 
               RowBox[{"Epilog", "\[Rule]", 
                RowBox[{"{", 
                 RowBox[{"defEpilog", ",", 
                  RowBox[{"(*", 
                   RowBox[{"Original", " ", "epilog"}], "*)"}], 
                  RowBox[{"(*", 
                   RowBox[{"Replot", " ", "button"}], "*)"}], 
                  RowBox[{"EventHandler", "[", 
                   RowBox[{
                    RowBox[{"Inset", "[", 
                    RowBox[{
                    RowBox[{"Graphics", "[", 
                    RowBox[{
                    RowBox[{"{", 
                    RowBox[{"Black", ",", 
                    RowBox[{"Dynamic", "@", 
                    RowBox[{"Opacity", "@", 
                    RowBox[{"If", "[", 
                    RowBox[{"overButton", ",", 
                    RowBox[{"If", "[", 
                    RowBox[{"pushedButton", ",", ".2", ",", ".1"}], "]"}], 
                    ",", "0"}], "]"}]}]}], ",", 
                    RowBox[{"Rectangle", "[", 
                    RowBox[{
                    RowBox[{"{", 
                    RowBox[{"0", ",", "0"}], "}"}], ",", 
                    RowBox[{"{", 
                    RowBox[{"1", ",", "1"}], "}"}], ",", 
                    RowBox[{"RoundingRadius", "\[Rule]", 
                    RowBox[{"Offset", "@", "5"}]}]}], "]"}], ",", "Black", 
                    ",", 
                    RowBox[{"Dynamic", "@", 
                    RowBox[{"Opacity", "@", 
                    RowBox[{"If", "[", 
                    RowBox[{"overButton", ",", 
                    RowBox[{"If", "[", 
                    RowBox[{"pushedButton", ",", ".9", ",", ".7"}], "]"}], 
                    ",", "0"}], "]"}]}]}], ",", 
                    RowBox[{"Text", "[", 
                    RowBox[{"\"\<Replot\>\"", ",", "Center"}], "]"}]}], "}"}],
                     ",", 
                    RowBox[{"ImageSize", "\[Rule]", "50"}], ",", 
                    RowBox[{"AspectRatio", "\[Rule]", 
                    RowBox[{"1", "/", "2"}]}], ",", 
                    RowBox[{"BaseStyle", "\[Rule]", 
                    RowBox[{"{", 
                    RowBox[{"FontFamily", "\[Rule]", "\"\<Helvetica\>\""}], 
                    "}"}]}]}], "]"}], ",", 
                    RowBox[{"Scaled", "@", 
                    RowBox[{"{", 
                    RowBox[{"1", ",", "1"}], "}"}]}], ",", 
                    RowBox[{"Scaled", "@", 
                    RowBox[{"{", 
                    RowBox[{"1", ",", "1"}], "}"}]}]}], "]"}], ",", 
                    RowBox[{"{", 
                    RowBox[{
                    RowBox[{"\"\<MouseEntered\>\"", "\[RuleDelayed]", 
                    RowBox[{"(", 
                    RowBox[{"overButton", "=", "True"}], ")"}]}], ",", 
                    RowBox[{"\"\<MouseExited\>\"", "\[RuleDelayed]", 
                    RowBox[{"(", 
                    RowBox[{"overButton", "=", "False"}], ")"}]}], ",", 
                    RowBox[{"\"\<MouseDown\>\"", "\[RuleDelayed]", 
                    RowBox[{"(", 
                    RowBox[{"pushedButton", "=", "True"}], ")"}]}], ",", 
                    RowBox[{"\"\<MouseUp\>\"", "\[RuleDelayed]", 
                    RowBox[{"(", 
                    RowBox[{
                    RowBox[{"pushedButton", "=", "False"}], ";", 
                    "\[IndentingNewLine]", 
                    RowBox[{"replot", "@", "plot"}]}], ")"}]}]}], "}"}]}], 
                   "]"}], ",", 
                  RowBox[{"(*", 
                   RowBox[{"Zoom", " ", "rectangle"}], "*)"}], 
                  RowBox[{"FaceForm", "@", 
                   RowBox[{"{", 
                    RowBox[{"Blue", ",", 
                    RowBox[{"Opacity", "@", ".1"}]}], "}"}]}], ",", 
                  RowBox[{"EdgeForm", "@", 
                   RowBox[{"{", 
                    RowBox[{"Thin", ",", "Dotted", ",", 
                    RowBox[{"Opacity", "@", ".5"}]}], "}"}]}], ",", 
                  RowBox[{"Dynamic", "@", "rect"}], ",", 
                  RowBox[{"(*", 
                   RowBox[{"Range", " ", "sliders"}], "*)"}], 
                  RowBox[{"Inset", "[", 
                   RowBox[{
                    RowBox[{"slider", "[", 
                    RowBox[{
                    RowBox[{"Dynamic", "@", "xPos"}], ",", 
                    RowBox[{"Dynamic", "@", "overHor"}], ",", 
                    RowBox[{"Dynamic", "@", "activeHor"}], ",", "size", ",", 
                    "True"}], "]"}], ",", 
                    RowBox[{"Scaled", "@", 
                    RowBox[{"{", 
                    RowBox[{".5", ",", "0"}], "}"}]}], ",", 
                    RowBox[{"Scaled", "@", 
                    RowBox[{"{", 
                    RowBox[{".5", ",", "0"}], "}"}]}]}], "]"}], ",", 
                  RowBox[{"Inset", "[", 
                   RowBox[{
                    RowBox[{"slider", "[", 
                    RowBox[{
                    RowBox[{"Dynamic", "@", "yPos"}], ",", 
                    RowBox[{"Dynamic", "@", "overVer"}], ",", 
                    RowBox[{"Dynamic", "@", "activeVer"}], ",", "size", ",", 
                    "False"}], "]"}], ",", 
                    RowBox[{"Scaled", "@", 
                    RowBox[{"{", 
                    RowBox[{"0", ",", ".5"}], "}"}]}], ",", 
                    RowBox[{"Scaled", "@", 
                    RowBox[{"{", 
                    RowBox[{"0", ",", ".5"}], "}"}]}]}], "]"}]}], "}"}]}]}], 
              "]"}], ",", 
             RowBox[{"(*", "cursor", "*)"}], 
             RowBox[{"Which", "[", 
              RowBox[{"dragged", ",", "\"\<FrameRisingResize\>\"", ",", 
               RowBox[{"CurrentValue", "@", "\"\<AltKey\>\""}], ",", 
               RowBox[{"Graphics", "[", 
                RowBox[{
                 RowBox[{"{", 
                  RowBox[{
                   RowBox[{"Antialiasing", "\[Rule]", "False"}], ",", 
                   RowBox[{"Line", "@", 
                    RowBox[{"{", 
                    RowBox[{
                    RowBox[{"{", 
                    RowBox[{
                    RowBox[{"-", "1"}], ",", "0"}], "}"}], ",", 
                    RowBox[{"{", 
                    RowBox[{"1", ",", "0"}], "}"}]}], "}"}]}], ",", 
                   RowBox[{"Line", "@", 
                    RowBox[{"{", 
                    RowBox[{
                    RowBox[{"{", 
                    RowBox[{"0", ",", 
                    RowBox[{"-", "1"}]}], "}"}], ",", 
                    RowBox[{"{", 
                    RowBox[{"0", ",", "1"}], "}"}]}], "}"}]}], ",", 
                   RowBox[{"Dynamic", "@", 
                    RowBox[{"Text", "[", 
                    RowBox[{
                    RowBox[{"Style", "[", 
                    RowBox[{"coord", ",", 
                    RowBox[{"FontFamily", "\[Rule]", "\"\<Helvetica\>\""}], 
                    ",", "8"}], "]"}], ",", 
                    RowBox[{"{", 
                    RowBox[{"0", ",", "0"}], "}"}], ",", 
                    RowBox[{"{", 
                    RowBox[{
                    RowBox[{"-", "1.1"}], ",", 
                    RowBox[{"-", "1"}]}], "}"}]}], "]"}]}]}], "}"}], ",", 
                 RowBox[{"ImageSize", "\[Rule]", "160"}], ",", 
                 RowBox[{"ImagePadding", "\[Rule]", "0"}], ",", 
                 RowBox[{"AspectRatio", "\[Rule]", "1"}], ",", 
                 RowBox[{"PlotRange", "\[Rule]", 
                  RowBox[{
                   RowBox[{"{", 
                    RowBox[{
                    RowBox[{"{", 
                    RowBox[{
                    RowBox[{"-", "1"}], ",", "1"}], "}"}], ",", 
                    RowBox[{"{", 
                    RowBox[{
                    RowBox[{"-", "1"}], ",", "1"}], "}"}]}], "}"}], "*", 
                   "8"}]}]}], "]"}], ",", 
               RowBox[{
                RowBox[{"CurrentValue", "@", "\"\<ControlKey\>\""}], "&&", 
                RowBox[{"!", "overHor"}], "&&", 
                RowBox[{"!", "overVer"}]}], ",", "\"\<ZoomView\>\"", ",", 
               RowBox[{
                RowBox[{"CurrentValue", "@", "\"\<ShiftKey\>\""}], "&&", 
                RowBox[{"!", "overHor"}], "&&", 
                RowBox[{"!", "overVer"}]}], ",", "\"\<DragGraphics\>\"", ",", 
               RowBox[{
                RowBox[{"CurrentValue", "@", "\"\<ShiftKey\>\""}], "&&", 
                RowBox[{"(", 
                 RowBox[{"overHor", "||", "overVer"}], ")"}]}], ",", 
               "\"\<LinkHand\>\"", ",", "True", ",", "Automatic"}], "]"}]}], 
            "]"}], ",", 
           RowBox[{"TrackedSymbols", "\[RuleDelayed]", 
            RowBox[{"{", 
             RowBox[{"dragged", ",", "plot"}], "}"}]}]}], "]"}], ",", 
         RowBox[{"{", 
          RowBox[{
           RowBox[{"\"\<MouseMoved\>\"", "\[RuleDelayed]", 
            RowBox[{"(", 
             RowBox[{"coord", "=", 
              RowBox[{"MousePosition", "@", "\"\<Graphics\>\""}]}], ")"}]}], 
           ",", 
           RowBox[{"\"\<MouseClicked\>\"", "\[RuleDelayed]", 
            RowBox[{"If", "[", 
             RowBox[{
              RowBox[{
               RowBox[{"!", "overButton"}], "&&", 
               RowBox[{
                RowBox[{"CurrentValue", "@", "\"\<MouseClickCount\>\""}], 
                "\[Equal]", "2"}]}], ",", 
              RowBox[{"reset", "[", "]"}]}], "]"}]}], ",", 
           RowBox[{"\"\<MouseDown\>\"", "\[RuleDelayed]", 
            RowBox[{"If", "[", 
             RowBox[{
              RowBox[{"(", 
               RowBox[{
                RowBox[{"!", "overHor"}], "||", 
                RowBox[{"!", "overVer"}]}], ")"}], ",", 
              RowBox[{
               RowBox[{"pt1", "=", 
                RowBox[{"MousePosition", "@", "\"\<Graphics\>\""}]}], ";", 
               "\[IndentingNewLine]", 
               RowBox[{"spt1", "=", 
                RowBox[{
                "MousePosition", "@", "\"\<GraphicsScaled\>\""}]}]}]}], 
             "]"}]}], ",", 
           RowBox[{"\"\<MouseDragged\>\"", "\[RuleDelayed]", 
            RowBox[{"If", "[", 
             RowBox[{
              RowBox[{
               RowBox[{"(", 
                RowBox[{
                 RowBox[{
                  RowBox[{"!", "overHor"}], "&&", 
                  RowBox[{"!", "overVer"}]}], "||", "dragged"}], ")"}], "&&", 
               RowBox[{"!", "overButton"}], "&&", 
               RowBox[{"!", "pushedButton"}], "&&", 
               RowBox[{"!", "activeVer"}], "&&", 
               RowBox[{"!", "activeHor"}]}], ",", 
              RowBox[{
               RowBox[{"dragged", "=", "True"}], ";", "\[IndentingNewLine]", 
               RowBox[{"Which", "[", 
                RowBox[{
                 RowBox[{"CurrentValue", "@", "\"\<ControlKey\>\""}], ",", 
                 RowBox[{
                  RowBox[{"dist", "=", 
                   RowBox[{
                    RowBox[{"Last", "@", 
                    RowBox[{
                    "MousePosition", "@", "\"\<GraphicsScaled\>\""}]}], "-", 
                    RowBox[{"Last", "@", "spt1"}]}]}], ";", 
                  "\[IndentingNewLine]", 
                  RowBox[{"range", "=", 
                   RowBox[{"zoomRange", "*", 
                    RowBox[{"10", "^", "dist"}]}]}]}], ",", 
                 RowBox[{"(*", 
                  RowBox[{
                  "\"\<GraphicsScaled\>\"", " ", "is", " ", "required", " ", 
                   "below", " ", "as", " ", "\"\<Graphics\>\"", " ", "would", 
                   " ", "give", " ", "jumpy", " ", "results", " ", "as", " ", 
                   "the", " ", "plot", " ", "range", " ", "is", " ", 
                   "constantly", " ", 
                   RowBox[{"modified", "."}]}], "*)"}], 
                 RowBox[{"CurrentValue", "@", "\"\<ShiftKey\>\""}], ",", 
                 RowBox[{
                  RowBox[{"pt2", "=", 
                   RowBox[{
                    RowBox[{"MapThread", "[", 
                    RowBox[{"Rescale", ",", 
                    RowBox[{"{", 
                    RowBox[{
                    RowBox[{"MousePosition", "@", "\"\<GraphicsScaled\>\""}], 
                    ",", 
                    RowBox[{"{", 
                    RowBox[{
                    RowBox[{"{", 
                    RowBox[{"0", ",", "1"}], "}"}], ",", 
                    RowBox[{"{", 
                    RowBox[{"0", ",", "1"}], "}"}]}], "}"}], ",", 
                    "zoomRange"}], "}"}]}], "]"}], "-", "pt1"}]}], ";", 
                  "\[IndentingNewLine]", 
                  RowBox[{"range", "=", 
                   RowBox[{"zoomRange", "-", "pt2"}]}]}], ",", "True", ",", 
                 RowBox[{
                  RowBox[{"pt2", "=", 
                   RowBox[{"MousePosition", "@", "\"\<Graphics\>\""}]}], ";", 
                  "\[IndentingNewLine]", 
                  RowBox[{"rect", "=", 
                   RowBox[{"If", "[", 
                    RowBox[{
                    RowBox[{"pt2", "===", "None"}], ",", 
                    RowBox[{"{", "}"}], ",", 
                    RowBox[{"Rectangle", "[", 
                    RowBox[{"pt1", ",", "pt2"}], "]"}]}], "]"}]}]}]}], 
                "]"}]}]}], "]"}]}], ",", 
           RowBox[{"\"\<MouseUp\>\"", "\[RuleDelayed]", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{"If", "[", 
               RowBox[{
                RowBox[{
                 RowBox[{"!", 
                  RowBox[{"CurrentValue", "@", "\"\<ShiftKey\>\""}]}], "&&", 
                 RowBox[{"!", 
                  RowBox[{"CurrentValue", "@", "\"\<ControlKey\>\""}]}], "&&",
                  "dragged", "&&", 
                 RowBox[{
                  RowBox[{"distance", "[", 
                   RowBox[{"pt1", ",", "pt2", ",", "zoomRange"}], "]"}], ">", 
                  "minDistance"}]}], ",", 
                RowBox[{
                 RowBox[{"rect", "=", 
                  RowBox[{"{", "}"}]}], ";", "\[IndentingNewLine]", 
                 RowBox[{"range", "=", 
                  RowBox[{"Transpose", "@", 
                   RowBox[{"Sort", "@", 
                    RowBox[{"{", 
                    RowBox[{"pt1", ",", "pt2"}], "}"}]}]}]}], ";"}]}], "]"}], 
              ";", 
              RowBox[{"zoomRange", "=", "range"}], ";", "\[IndentingNewLine]", 
              RowBox[{"dragged", "=", "False"}]}], ")"}]}]}], "}"}], ",", 
         RowBox[{"PassEventsDown", "\[Rule]", "True"}], ",", 
         RowBox[{"PassEventsUp", "\[Rule]", "True"}]}], "]"}]}]}], ",", 
     RowBox[{"Initialization", "\[RuleDelayed]", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{
         RowBox[{"slider", "[", 
          RowBox[{
           RowBox[{"Dynamic", "[", "pos_", "]"}], ",", 
           RowBox[{"Dynamic", "[", "over_", "]"}], ",", 
           RowBox[{"Dynamic", "[", "active_", "]"}], ",", "size_", ",", 
           "hor_"}], "]"}], ":=", 
         RowBox[{"EventHandler", "[", 
          RowBox[{
           RowBox[{"Dynamic", "[", 
            RowBox[{
             RowBox[{"LocatorPane", "[", 
              RowBox[{
               RowBox[{"If", "[", 
                RowBox[{"hor", ",", 
                 RowBox[{"Dynamic", "[", 
                  RowBox[{
                   RowBox[{"{", 
                    RowBox[{"pos", ",", ".5"}], "}"}], ",", 
                   RowBox[{
                    RowBox[{"(", 
                    RowBox[{
                    RowBox[{"pos", "=", 
                    RowBox[{"Clip", "[", 
                    RowBox[{
                    RowBox[{"First", "@", "#"}], ",", 
                    RowBox[{"{", 
                    RowBox[{"0", ",", "1"}], "}"}]}], "]"}]}], ";", 
                    "\[IndentingNewLine]", 
                    RowBox[{
                    RowBox[{"range", "[", 
                    RowBox[{"[", "1", "]"}], "]"}], "=", 
                    RowBox[{"If", "[", 
                    RowBox[{
                    RowBox[{"CurrentValue", "@", "\"\<ShiftKey\>\""}], ",", 
                    RowBox[{
                    RowBox[{"First", "@", "zoomRange"}], "+", 
                    RowBox[{
                    RowBox[{"(", 
                    RowBox[{"pos", "-", ".5"}], ")"}], "*", 
                    RowBox[{"Abs", "[", 
                    RowBox[{"Subtract", "@@", 
                    RowBox[{"First", "@", "zoomRange"}]}], "]"}], "*", 
                    "panFactor"}]}], ",", 
                    RowBox[{
                    RowBox[{"First", "@", "zoomRange"}], "*", 
                    RowBox[{"10", "^", 
                    RowBox[{"(", 
                    RowBox[{
                    RowBox[{"(", 
                    RowBox[{"pos", "-", ".5"}], ")"}], "*", "2"}], 
                    ")"}]}]}]}], "]"}]}]}], ")"}], "&"}]}], "]"}], ",", 
                 RowBox[{"Dynamic", "[", 
                  RowBox[{
                   RowBox[{"{", 
                    RowBox[{".5", ",", "pos"}], "}"}], ",", 
                   RowBox[{
                    RowBox[{"(", 
                    RowBox[{
                    RowBox[{"pos", "=", 
                    RowBox[{"Clip", "[", 
                    RowBox[{
                    RowBox[{"Last", "@", "#"}], ",", 
                    RowBox[{"{", 
                    RowBox[{"0", ",", "1"}], "}"}]}], "]"}]}], ";", 
                    "\[IndentingNewLine]", 
                    RowBox[{
                    RowBox[{"range", "[", 
                    RowBox[{"[", "2", "]"}], "]"}], "=", 
                    RowBox[{"If", "[", 
                    RowBox[{
                    RowBox[{"CurrentValue", "@", "\"\<ShiftKey\>\""}], ",", 
                    RowBox[{
                    RowBox[{"Last", "@", "zoomRange"}], "+", 
                    RowBox[{
                    RowBox[{"(", 
                    RowBox[{"pos", "-", ".5"}], ")"}], "*", 
                    RowBox[{"Abs", "[", 
                    RowBox[{"Subtract", "@@", 
                    RowBox[{"Last", "@", "zoomRange"}]}], "]"}], "*", 
                    "panFactor"}]}], ",", 
                    RowBox[{
                    RowBox[{"Last", "@", "zoomRange"}], "*", 
                    RowBox[{"10", "^", 
                    RowBox[{"(", 
                    RowBox[{
                    RowBox[{"(", 
                    RowBox[{"pos", "-", ".5"}], ")"}], "*", "2"}], 
                    ")"}]}]}]}], "]"}]}]}], ")"}], "&"}]}], "]"}]}], "]"}], 
               ",", 
               RowBox[{"(*", 
                RowBox[{
                "Locator", " ", "below", " ", "is", " ", "required", " ", 
                 "inside", " ", "Graphics", " ", "as", " ", "LocatorPane", 
                 " ", "inside", " ", "Inset", " ", "looses", " ", "the", " ", 
                 "crosshair", " ", 
                 RowBox[{"image", "."}]}], "*)"}], 
               RowBox[{"Graphics", "[", 
                RowBox[{
                 RowBox[{"{", 
                  RowBox[{"Black", ",", 
                   RowBox[{"Opacity", "@", 
                    RowBox[{"If", "[", 
                    RowBox[{
                    RowBox[{
                    RowBox[{"(", 
                    RowBox[{"over", "||", "active"}], ")"}], "&&", 
                    RowBox[{"!", "dragged"}]}], ",", ".1", ",", "0"}], 
                    "]"}]}], ",", 
                   RowBox[{"Rectangle", "[", 
                    RowBox[{
                    RowBox[{"{", 
                    RowBox[{"0", ",", "0"}], "}"}], ",", 
                    RowBox[{"{", 
                    RowBox[{"1", ",", "1"}], "}"}], ",", 
                    RowBox[{"RoundingRadius", "\[Rule]", 
                    RowBox[{"Offset", "@", "5"}]}]}], "]"}], ",", 
                   RowBox[{"Opacity", "@", "1"}], ",", 
                   RowBox[{"Locator", "[", 
                    RowBox[{
                    RowBox[{"Dynamic", "@", 
                    RowBox[{"If", "[", 
                    RowBox[{"hor", ",", 
                    RowBox[{"{", 
                    RowBox[{"pos", ",", ".5"}], "}"}], ",", 
                    RowBox[{"{", 
                    RowBox[{".5", ",", "pos"}], "}"}]}], "]"}]}], ",", 
                    RowBox[{"Appearance", "\[Rule]", 
                    RowBox[{"If", "[", 
                    RowBox[{
                    RowBox[{
                    RowBox[{"(", 
                    RowBox[{"over", "||", "active"}], ")"}], "&&", 
                    RowBox[{"!", "dragged"}]}], ",", "Automatic", ",", 
                    "None"}], "]"}]}]}], "]"}]}], "}"}], ",", 
                 RowBox[{"ImagePadding", "\[Rule]", "0"}], ",", 
                 RowBox[{"PlotRangePadding", "\[Rule]", "0"}], ",", 
                 RowBox[{"ImageSize", "\[Rule]", 
                  RowBox[{
                   RowBox[{"If", "[", 
                    RowBox[{"hor", ",", 
                    RowBox[{
                    RowBox[{
                    RowBox[{"#", "*", 
                    RowBox[{"{", 
                    RowBox[{"1", ",", "0"}], "}"}]}], "+", 
                    RowBox[{"{", 
                    RowBox[{"0", ",", "20"}], "}"}]}], "&"}], ",", 
                    RowBox[{
                    RowBox[{
                    RowBox[{"#", "*", 
                    RowBox[{"{", 
                    RowBox[{"0", ",", "1"}], "}"}]}], "+", 
                    RowBox[{"{", 
                    RowBox[{"20", ",", "0"}], "}"}]}], "&"}]}], "]"}], "@", 
                   "size"}]}], ",", 
                 RowBox[{"AspectRatio", "\[Rule]", 
                  RowBox[{"If", "[", 
                   RowBox[{"hor", ",", 
                    RowBox[{
                    RowBox[{"1", "/", "sliderRatio"}], "/", 
                    RowBox[{"First", "@", "size"}]}], ",", 
                    RowBox[{"sliderRatio", "*", 
                    RowBox[{"Last", "@", "size"}]}]}], "]"}]}]}], "]"}], ",", 
               RowBox[{"Enabled", "\[Rule]", 
                RowBox[{"(", 
                 RowBox[{"over", "||", "active"}], ")"}]}], ",", 
               RowBox[{"Appearance", "\[Rule]", "None"}], ",", 
               RowBox[{"FrameMargins", "\[Rule]", "0"}], ",", 
               RowBox[{"ContentPadding", "\[Rule]", "False"}]}], "]"}], ",", 
             RowBox[{"TrackedSymbols", "\[RuleDelayed]", 
              RowBox[{"{", 
               RowBox[{"pos", ",", "over", ",", "active", ",", "dragged"}], 
               "}"}]}]}], "]"}], ",", 
           RowBox[{"{", 
            RowBox[{
             RowBox[{"\"\<MouseDown\>\"", "\[RuleDelayed]", 
              RowBox[{"(", 
               RowBox[{"active", "=", "True"}], ")"}]}], ",", 
             RowBox[{"\"\<MouseUp\>\"", "\[RuleDelayed]", 
              RowBox[{"(", 
               RowBox[{"active", "=", "False"}], ")"}]}], ",", 
             RowBox[{"\"\<MouseEntered\>\"", "\[RuleDelayed]", 
              RowBox[{"(", 
               RowBox[{"over", "=", 
                RowBox[{"!", 
                 RowBox[{"If", "[", 
                  RowBox[{"hor", ",", "activeVer", ",", "activeHor"}], 
                  "]"}]}]}], ")"}]}], ",", 
             RowBox[{"(*", 
              RowBox[{
               RowBox[{
               "The", " ", "value", " ", "of", " ", "`active`", " ", "must", 
                " ", "be", " ", "public", " ", "so", " ", "that", " ", "one", 
                " ", "slider", " ", "can", " ", "check", " ", "whether", " ", 
                "the", " ", "other", " ", "is", " ", "dragged"}], ",", 
               RowBox[{
               "not", " ", "to", " ", "pop", " ", "up", " ", "if", " ", "the",
                 " ", "other", " ", "locator", " ", "is", " ", "moved", " ", 
                "over", " ", "this", " ", 
                RowBox[{"slider", "'"}], "s", " ", 
                RowBox[{"area", "."}]}]}], "*)"}], 
             RowBox[{"\"\<MouseExited\>\"", "\[RuleDelayed]", 
              RowBox[{"(", 
               RowBox[{"over", "=", "False"}], ")"}]}]}], "}"}], ",", 
           RowBox[{"PassEventsDown", "\[Rule]", "True"}]}], "]"}]}], ";", 
        "\[IndentingNewLine]", 
        RowBox[{
         RowBox[{"replot", "[", "p_", "]"}], ":=", 
         RowBox[{"Module", "[", 
          RowBox[{
           RowBox[{"{", "temp", "}"}], ",", 
           RowBox[{"Switch", "[", 
            RowBox[{"head", ",", "Plot", ",", 
             RowBox[{
              RowBox[{"temp", "=", 
               RowBox[{"ReplacePart", "[", 
                RowBox[{"held", ",", 
                 RowBox[{"{", 
                  RowBox[{
                   RowBox[{"{", 
                    RowBox[{"1", ",", "2"}], "}"}], "\[Rule]", 
                   RowBox[{"{", 
                    RowBox[{
                    RowBox[{"held", "[", 
                    RowBox[{"[", 
                    RowBox[{"1", ",", "2", ",", "1"}], "]"}], "]"}], ",", 
                    RowBox[{"range", "[", 
                    RowBox[{"[", 
                    RowBox[{"1", ",", "1"}], "]"}], "]"}], ",", 
                    RowBox[{"range", "[", 
                    RowBox[{"[", 
                    RowBox[{"1", ",", "2"}], "]"}], "]"}]}], "}"}]}], "}"}]}],
                 "]"}]}], ";", "\[IndentingNewLine]", 
              RowBox[{"(*", 
               RowBox[{
               "Replace", " ", "plot", " ", "range", " ", "or", " ", "insert",
                 " ", "if", " ", "nonexistent"}], "*)"}], 
              RowBox[{"plot", "=", 
               RowBox[{"ReleaseHold", "@", 
                RowBox[{"If", "[", 
                 RowBox[{
                  RowBox[{"MemberQ", "[", 
                   RowBox[{"temp", ",", "PlotRange", ",", "Infinity"}], "]"}],
                   ",", 
                  RowBox[{"temp", "/.", 
                   RowBox[{"{", 
                    RowBox[{
                    RowBox[{"_", "[", 
                    RowBox[{"PlotRange", ",", "_"}], "]"}], "\[Rule]", 
                    RowBox[{"(", 
                    RowBox[{"PlotRange", "\[Rule]", "range"}], ")"}]}], 
                    "}"}]}], ",", 
                  RowBox[{"Insert", "[", 
                   RowBox[{"temp", ",", 
                    RowBox[{"PlotRange", "\[Rule]", "range"}], ",", 
                    RowBox[{"{", 
                    RowBox[{"1", ",", 
                    RowBox[{"-", "1"}]}], "}"}]}], "]"}]}], "]"}]}]}]}], ",", 
             RowBox[{
             "DensityPlot", "|", "ContourPlot", "|", "RegionPlot", "|", 
              "StreamPlot", "|", "StreamDensityPlot", "|", "VectorPlot", "|", 
              "VectorDensityPlot", "|", "LineIntegralConvolutionPlot"}], ",", 
             RowBox[{
              RowBox[{"temp", "=", 
               RowBox[{"ReplacePart", "[", 
                RowBox[{"held", ",", 
                 RowBox[{"{", 
                  RowBox[{
                   RowBox[{
                    RowBox[{"{", 
                    RowBox[{"1", ",", "2"}], "}"}], "\[Rule]", 
                    RowBox[{"{", 
                    RowBox[{
                    RowBox[{"held", "[", 
                    RowBox[{"[", 
                    RowBox[{"1", ",", "2", ",", "1"}], "]"}], "]"}], ",", 
                    RowBox[{"range", "[", 
                    RowBox[{"[", 
                    RowBox[{"1", ",", "1"}], "]"}], "]"}], ",", 
                    RowBox[{"range", "[", 
                    RowBox[{"[", 
                    RowBox[{"1", ",", "2"}], "]"}], "]"}]}], "}"}]}], ",", 
                   RowBox[{
                    RowBox[{"{", 
                    RowBox[{"1", ",", "3"}], "}"}], "\[Rule]", 
                    RowBox[{"{", 
                    RowBox[{
                    RowBox[{"held", "[", 
                    RowBox[{"[", 
                    RowBox[{"1", ",", "3", ",", "1"}], "]"}], "]"}], ",", 
                    RowBox[{"range", "[", 
                    RowBox[{"[", 
                    RowBox[{"2", ",", "1"}], "]"}], "]"}], ",", 
                    RowBox[{"range", "[", 
                    RowBox[{"[", 
                    RowBox[{"2", ",", "2"}], "]"}], "]"}]}], "}"}]}]}], 
                  "}"}]}], "]"}]}], ";", "\[IndentingNewLine]", 
              RowBox[{"plot", "=", 
               RowBox[{"ReleaseHold", "@", 
                RowBox[{"If", "[", 
                 RowBox[{
                  RowBox[{"MemberQ", "[", 
                   RowBox[{"temp", ",", "PlotRange", ",", "Infinity"}], "]"}],
                   ",", 
                  RowBox[{"temp", "/.", 
                   RowBox[{"{", 
                    RowBox[{
                    RowBox[{"_", "[", 
                    RowBox[{"PlotRange", ",", "_"}], "]"}], "\[Rule]", 
                    RowBox[{"(", 
                    RowBox[{"PlotRange", "\[Rule]", "range"}], ")"}]}], 
                    "}"}]}], ",", 
                  RowBox[{"Insert", "[", 
                   RowBox[{"temp", ",", 
                    RowBox[{"PlotRange", "\[Rule]", "range"}], ",", 
                    RowBox[{"{", 
                    RowBox[{"1", ",", 
                    RowBox[{"-", "1"}]}], "}"}]}], "]"}]}], "]"}]}]}], ";"}], 
             ",", "_", ",", "Null"}], "]"}]}], "]"}]}], ";", 
        "\[IndentingNewLine]", 
        RowBox[{
         RowBox[{"reset", "[", "]"}], ":=", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"range", "=", 
            RowBox[{"zoomRange", "=", "defRange"}]}], ";", 
           RowBox[{"xPos", "=", 
            RowBox[{"yPos", "=", ".5"}]}]}], ")"}]}], ";", 
        "\[IndentingNewLine]", 
        RowBox[{"(*", 
         RowBox[{
          RowBox[{
          "Quiet", " ", "suppresses", " ", "error", " ", "message", " ", 
           "when", " ", 
           RowBox[{"Rescale", "[", 
            RowBox[{"y", ",", 
             RowBox[{"{", 
              RowBox[{"x", ",", "x"}], "}"}]}], "]"}], " ", "where", " ", 
           "y"}], "\[NotEqual]", 
          RowBox[{
           RowBox[{"x", ".", "EuclideanDistance"}], " ", "then", " ", 
           "defaults", " ", "to", " ", "Infinity", " ", "which", " ", "is", 
           " ", 
           RowBox[{"fine", "."}]}]}], "*)"}], 
        RowBox[{
         RowBox[{"distance", "[", 
          RowBox[{"p1_", ",", "p2_", ",", 
           RowBox[{"{", 
            RowBox[{"r1_", ",", "r2_"}], "}"}]}], "]"}], ":=", 
         RowBox[{"Quiet", "[", 
          RowBox[{"N", "@", 
           RowBox[{"EuclideanDistance", "[", 
            RowBox[{
             RowBox[{"Rescale", "[", 
              RowBox[{"p1", ",", "r1"}], "]"}], ",", 
             RowBox[{"Rescale", "[", 
              RowBox[{"p2", ",", "r2"}], "]"}]}], "]"}]}], "]"}]}], ";"}], 
       ")"}]}]}], "]"}]}], ";"}], "\n"}], "Input",
 CellChangeTimes->{{3.6902655649832325`*^9, 3.690265578087485*^9}, 
   3.6902674480999365`*^9, 3.6902675212394743`*^9, {3.6902682214555125`*^9, 
   3.690268221675517*^9}, {3.6902682563842325`*^9, 3.690268281554736*^9}, {
   3.690268313543395*^9, 3.6902683166234565`*^9}, {3.690268348654097*^9, 
   3.6902684147454453`*^9}}]
}, Open  ]]
}, Closed]],

Cell[CellGroupData[{

Cell["\<\
Get survey files:
psb.survey  <-- Output from 2016 
PSB_BEAM_COORDINATES.xlsx\
\>", "Title",
 CellChangeTimes->{{3.496374734535734*^9, 3.4963747405356956`*^9}, {
   3.4963936361647625`*^9, 3.496393636321012*^9}, {3.496393666961441*^9, 
   3.496393679164488*^9}, {3.496470149404563*^9, 3.496470183076223*^9}, {
   3.520161272116864*^9, 3.520161286365679*^9}, 3.6053439373064137`*^9, {
   3.605442764608918*^9, 3.6054427655761366`*^9}, {3.605500595662073*^9, 
   3.6055005958336763`*^9}, {3.6056367730729856`*^9, 
   3.6056367743391128`*^9}, {3.606207166252028*^9, 3.606207166829239*^9}, {
   3.606207732767044*^9, 3.606207733375456*^9}, {3.6089803778063817`*^9, 
   3.6089804104890103`*^9}, 3.609137412887022*^9, {3.60914709067869*^9, 
   3.609147152395059*^9}, {3.6862871980554943`*^9, 3.686287198775509*^9}, {
   3.6862892371712065`*^9, 3.686289241541294*^9}, {3.690263397641119*^9, 
   3.690263460463527*^9}},
 Background->RGBColor[1., 0., 1.]],

Cell[BoxData[{
 RowBox[{
  RowBox[{"fileDB", "=", 
   RowBox[{
    RowBox[{"Import", "[", 
     RowBox[{
     "\"\<http://project-ps-optics.web.cern.ch/project-PS-optics/cps/Psb/post_\
EYETS_2016_2017/cmd/PSB_BEAM_COORDINATES_08DEC2016.xlsx\>\"", ",", 
      "\"\<Data\>\""}], "]"}], "[", 
    RowBox[{"[", "1", "]"}], "]"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"fileMADX", "=", 
   RowBox[{"Import", "[", 
    RowBox[{
    "\"\<http://project-ps-optics.web.cern.ch/project-PS-optics/cps/Psb/post_\
EYETS_2016_2017/output/psb.survey\>\"", ",", "\"\<Table\>\""}], "]"}]}], 
  ";"}], "\[IndentingNewLine]"}], "Input",
 CellChangeTimes->{{3.4961631541260786`*^9, 3.496163220906489*^9}, {
   3.4961632949469433`*^9, 3.4961633115290456`*^9}, {3.496163389683527*^9, 
   3.4961634150436826`*^9}, {3.4961637622778444`*^9, 3.4961637639578547`*^9}, 
   3.4964701125922985`*^9, 3.520160851854958*^9, {3.6053439476026115`*^9, 
   3.6053439482266235`*^9}, {3.6054427733762865`*^9, 
   3.6054427737350936`*^9}, {3.605463767118492*^9, 3.6054637672900953`*^9}, {
   3.6054638288020782`*^9, 3.6054638289892817`*^9}, 3.605499412504914*^9, 
   3.6055006002329607`*^9, {3.605636671811861*^9, 3.605636675091189*^9}, {
   3.6056367081204915`*^9, 3.605636716157295*^9}, {3.605637183797571*^9, 
   3.6056371876664457`*^9}, {3.6062071600275087`*^9, 3.606207160979127*^9}, {
   3.6062076640945234`*^9, 3.6062077164647307`*^9}, {3.6089804263701153`*^9, 
   3.608980427056529*^9}, {3.6091373348236217`*^9, 3.6091373350576243`*^9}, 
   3.6091374129182224`*^9, {3.609137449781495*^9, 3.6091374499530973`*^9}, {
   3.609140925495953*^9, 3.609140933358656*^9}, {3.609150242230831*^9, 
   3.6091502651635666`*^9}, {3.609150300139887*^9, 3.609150335911834*^9}, {
   3.609747634460168*^9, 3.609747647517619*^9}, {3.6097477581393466`*^9, 
   3.6097477592001667`*^9}, {3.6097479418796797`*^9, 
   3.6097479633612933`*^9}, {3.6097501096977215`*^9, 3.609750110727341*^9}, {
   3.60975050524693*^9, 3.609750566119301*^9}, {3.68628720783569*^9, 
   3.6862872185659046`*^9}, {3.6862893009844847`*^9, 3.68628932624499*^9}, {
   3.6862893639797487`*^9, 3.6862893818201056`*^9}, 3.686289414474771*^9, {
   3.6862894556955957`*^9, 3.6862895616997724`*^9}, {3.6902629602399077`*^9, 
   3.6902629960425963`*^9}, {3.6902630370245843`*^9, 
   3.6902630711892414`*^9}, {3.6902631679735026`*^9, 3.690263170063943*^9}, {
   3.6902632922142916`*^9, 3.6902632971127863`*^9}, {3.6902640879718056`*^9, 
   3.69026412720656*^9}, {3.6902641933362317`*^9, 3.6902642066900883`*^9}, {
   3.690264314971771*^9, 3.690264315814187*^9}}],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"fileDB", "[", 
  RowBox[{"[", 
   RowBox[{"Range", "[", 
    RowBox[{"2", ",", "4"}], "]"}], "]"}], "]"}]], "Input",
 CellChangeTimes->{{3.6902630812202344`*^9, 3.690263125259881*^9}, {
  3.690263189018307*^9, 3.690263189361514*^9}, {3.6902632302811007`*^9, 
  3.6902632506394925`*^9}, {3.6902633021828833`*^9, 3.6902633503098087`*^9}, {
  3.6902640971603823`*^9, 3.690264097519189*^9}}],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{"\<\"BR.BHZ.11.E\"\>", ",", "1879.49884`", ",", "2110.85661`", 
     ",", "2433.66`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"\<\"BR.BHZ.11.S\"\>", ",", "1878.74378`", ",", "2112.28435`", 
     ",", "2433.66`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"\<\"BR.QFO.11.E\"\>", ",", "1878.52847`", ",", "2112.60968`", 
     ",", "2433.66`"}], "}"}]}], "}"}]], "Output",
 CellChangeTimes->{{3.6902630954633083`*^9, 3.69026312574349*^9}, {
   3.6902631821541753`*^9, 3.690263189642319*^9}, {3.6902632332139573`*^9, 
   3.6902632511231017`*^9}, {3.6902633037741137`*^9, 3.6902633511210246`*^9}, 
   3.690265724839507*^9}]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"fileMADX", "[", 
  RowBox[{"[", 
   RowBox[{"Range", "[", 
    RowBox[{"1", ",", "9"}], "]"}], "]"}], "]"}]], "Input",
 CellChangeTimes->{{3.690264219310731*^9, 3.690264282772752*^9}, {
  3.690264320806283*^9, 3.6902643436295223`*^9}}],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{"\<\"@\"\>", ",", "\<\"NAME\"\>", ",", "\<\"%06s\"\>", 
     ",", "\<\"SURVEY\"\>"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"\<\"@\"\>", ",", "\<\"TYPE\"\>", ",", "\<\"%06s\"\>", 
     ",", "\<\"SURVEY\"\>"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"\<\"@\"\>", ",", "\<\"TITLE\"\>", ",", "\<\"%15s\"\>", 
     ",", "\<\"BOOSTER lattice\"\>"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"\<\"@\"\>", ",", "\<\"ORIGIN\"\>", ",", "\<\"%16s\"\>", 
     ",", "\<\"5.02.12 Linux 64\"\>"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"\<\"@\"\>", ",", "\<\"DATE\"\>", ",", "\<\"%08s\"\>", 
     ",", "\<\"08/12/16\"\>"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"\<\"@\"\>", ",", "\<\"TIME\"\>", ",", "\<\"%08s\"\>", 
     ",", "\<\"14.53.59\"\>"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"\<\"*\"\>", ",", "\<\"NAME\"\>", ",", "\<\"S\"\>", 
     ",", "\<\"L\"\>", ",", "\<\"Z\"\>", ",", "\<\"X\"\>", ",", "\<\"Y\"\>"}],
     "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"\<\"$\"\>", ",", "\<\"%s\"\>", ",", "\<\"%le\"\>", 
     ",", "\<\"%le\"\>", ",", "\<\"%le\"\>", ",", "\<\"%le\"\>", 
     ",", "\<\"%le\"\>"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"\<\"PSB4$START\"\>", ",", "0.`", ",", "0.`", ",", 
     "1880.5034756326`", ",", "2108.4001001986`", ",", "2433.66`"}], "}"}]}], 
  "}"}]], "Output",
 CellChangeTimes->{{3.6902642250360413`*^9, 3.6902642836619687`*^9}, {
   3.6902643283568287`*^9, 3.6902643445655403`*^9}, 3.6902657248707075`*^9}]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Length", "[", "fileDB", "]"}]], "Input",
 CellChangeTimes->{{3.6091374066157417`*^9, 3.609137441419788*^9}, {
  3.690263360028796*^9, 3.690263362946052*^9}, {3.690264100748451*^9, 
  3.6902641010760574`*^9}}],

Cell[BoxData["271"], "Output",
 CellChangeTimes->{
  3.421509079424489*^9, 3.4215098397661853`*^9, 3.421550840081599*^9, 
   3.4215565460279284`*^9, 3.4222525660355563`*^9, 3.422253048435443*^9, 
   3.4222533198781548`*^9, 3.4222546508412085`*^9, {3.4222549132681227`*^9, 
   3.4222549337060146`*^9}, 3.422255283306477*^9, 3.4222554079807453`*^9, 
   3.4222556135003166`*^9, 3.422265882978737*^9, 3.4222864853754845`*^9, 
   3.4223314955158286`*^9, 3.4223330096110287`*^9, 3.4223437355249166`*^9, 
   3.422547757140625*^9, 3.4226079237704754`*^9, 3.422690866101122*^9, 
   3.422709430645871*^9, 3.422724790109375*^9, 3.4227663627051735`*^9, 
   3.423216930520055*^9, 3.4232231864800167`*^9, 3.4232233249478803`*^9, 
   3.423223428369093*^9, 3.4232234645876117`*^9, 3.4232236579613743`*^9, 
   3.424600473964882*^9, 3.4246011880922117`*^9, 3.4246022460505257`*^9, 
   3.424602307402046*^9, 3.424668967200333*^9, 3.424748162409685*^9, 
   3.425791939509286*^9, 3.4258140667929583`*^9, 3.4258146046404195`*^9, 
   3.4258147854905343`*^9, 3.425874095967534*^9, 3.425877427791161*^9, 
   3.425905680653373*^9, 3.4259668677401233`*^9, 3.4262447789761515`*^9, 
   3.4262471964867387`*^9, 3.426259603953125*^9, 3.42631467691314*^9, 
   3.4961638029180937`*^9, 3.496199738927723*^9, 3.4963265584569006`*^9, 
   3.4963742684137173`*^9, 3.496375879200283*^9, 3.496392533609319*^9, 
   3.496392640514885*^9, 3.4963928959820004`*^9, 3.496393110215004*^9, 
   3.4963960653523407`*^9, 3.4963962818665805`*^9, 3.49646993424969*^9, 
   3.4964704699650116`*^9, 3.4964709782586336`*^9, 3.4964711163827496`*^9, 
   3.4964712409757023`*^9, 3.4964722957502017`*^9, 3.496472541608003*^9, 
   3.496473019386195*^9, 3.4964730684171314`*^9, 3.496474369389532*^9, 
   3.496474429719494*^9, 3.4964745142602687`*^9, 3.496474935257598*^9, 
   3.4964762789845643`*^9, 3.5201608580888147`*^9, 3.520161504160057*^9, 
   3.605344022967661*^9, 3.6053548019351482`*^9, 3.605374672909955*^9, 
   3.605438058920065*^9, 3.6054389357057395`*^9, 3.6054391940467076`*^9, 
   3.6054427790391955`*^9, 3.6054592121407275`*^9, 3.6054622046921167`*^9, 
   3.6054637759794626`*^9, 3.6054638354010053`*^9, 3.6054988600202894`*^9, 
   3.605499417387808*^9, 3.6055006078615074`*^9, 3.605636859251603*^9, 
   3.605637195638199*^9, 3.6057167407877636`*^9, 3.6062071813999195`*^9, 
   3.6062077626884193`*^9, 3.606208071418357*^9, 3.6065545544372053`*^9, 
   3.606555457106805*^9, 3.6065566572535267`*^9, 3.606556797908971*^9, 
   3.606557550813583*^9, 3.6065579016863327`*^9, 3.606559383372429*^9, 
   3.6065594515769405`*^9, 3.606561758228883*^9, 3.6066293418138433`*^9, {
   3.6066323125663185`*^9, 3.606632317074748*^9}, 3.606633797471059*^9, 
   3.606637666330842*^9, 3.6066376974179945`*^9, 3.6066377494796104`*^9, 
   3.6066378925973663`*^9, 3.6066380009065437`*^9, 3.6066380601356792`*^9, 
   3.606638119922471*^9, 3.606638179642716*^9, 3.606638681981094*^9, 
   3.606638775988055*^9, 3.6066442637643337`*^9, 3.6066446746949744`*^9, 
   3.6089804932026043`*^9, {3.609137432309271*^9, 3.6091374610293727`*^9}, 
   3.609150425022419*^9, 3.6091507580179663`*^9, 3.609153256154098*^9, 
   3.6091535982846303`*^9, 3.6097464355607104`*^9, 3.6097502113180757`*^9, 
   3.686287268860914*^9, 3.686455755017836*^9, 3.6864561869260626`*^9, 
   3.6902633639444714`*^9, 3.690265724901908*^9}]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Length", "[", "fileMADX", "]"}]], "Input",
 CellChangeTimes->{{3.6902643677007847`*^9, 3.6902643741749096`*^9}}],

Cell[BoxData["520"], "Output",
 CellChangeTimes->{3.690264374627318*^9, 3.690265724933109*^9}]
}, Open  ]],

Cell[BoxData[{
 RowBox[{
  RowBox[{"XYDB", "  ", "=", 
   RowBox[{"GetXYDB", "[", 
    RowBox[{"fileDB", "[", 
     RowBox[{"[", 
      RowBox[{"Range", "[", 
       RowBox[{"2", ",", "271"}], "]"}], "]"}], "]"}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"XYMADX", "=", 
   RowBox[{"GetXYMADX", "[", 
    RowBox[{"fileMADX", "[", 
     RowBox[{"[", 
      RowBox[{"Range", "[", 
       RowBox[{"9", ",", "520"}], "]"}], "]"}], "]"}], "]"}]}], 
  ";"}]}], "Input",
 CellChangeTimes->{
  3.605343985948149*^9, 3.6053440208616204`*^9, {3.605344761373888*^9, 
   3.605344819984215*^9}, {3.60534550489208*^9, 3.6053455157966895`*^9}, {
   3.6053472999670105`*^9, 3.6053473015894423`*^9}, 3.605498879707868*^9, {
   3.605636918308508*^9, 3.605636927952651*^9}, 3.605637137090273*^9, 
   3.6056374053842325`*^9, {3.6056378239870825`*^9, 3.605637845562297*^9}, {
   3.6065554668101916`*^9, 3.6065554677930107`*^9}, {3.608980549940895*^9, 
   3.6089805509705153`*^9}, {3.608980697036124*^9, 3.6089806975041327`*^9}, 
   3.60913741274662*^9, {3.609141259520854*^9, 3.609141344869548*^9}, {
   3.6091504503577933`*^9, 3.609150465942793*^9}, {3.609152166661189*^9, 
   3.609152168049651*^9}, {3.60915224357264*^9, 3.609152248206048*^9}, {
   3.6862873776430893`*^9, 3.686287378193101*^9}, {3.686287442098383*^9, 
   3.686287442788397*^9}, {3.690263676212686*^9, 3.6902636856664677`*^9}, {
   3.690263895963313*^9, 3.690263905370294*^9}, {3.690264036350413*^9, 
   3.6902640660689845`*^9}, {3.6902644203985987`*^9, 
   3.6902644569189005`*^9}, {3.6902657442930813`*^9, 3.690265749191575*^9}}],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"PlotExplorer", "@", 
  RowBox[{"ListLinePlot", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{"XYDB", ",", "XYMADX"}], "}"}], ",", 
    RowBox[{"PlotStyle", "\[Rule]", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"{", "Blue", "}"}], ",", 
       RowBox[{"{", 
        RowBox[{"Red", ",", "Dashed"}], "}"}]}], "}"}]}], ",", 
    RowBox[{"PlotLegends", "\[Rule]", 
     RowBox[{"{", 
      RowBox[{"\"\<Survey Database\>\"", ",", " ", "\"\<Survey MADX\>\""}], 
      "}"}]}], ",", " ", 
    RowBox[{"PlotRange", "\[Rule]", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"{", 
        RowBox[{"1831", ",", "1882"}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{"2075", ",", "2126"}], "}"}]}], "}"}]}], ",", 
    RowBox[{"ImageSize", "\[Rule]", "1000"}], ",", 
    RowBox[{"AspectRatio", "\[Rule]", "1"}]}], "]"}]}]], "Input",
 CellChangeTimes->{{3.690265591768948*^9, 3.6902655984146757`*^9}, {
  3.6902659188448377`*^9, 3.6902659391720285`*^9}, {3.690266026861315*^9, 
  3.6902660944560533`*^9}, {3.690268074898522*^9, 3.6902681372978067`*^9}, {
  3.690268181972721*^9, 3.6902681847627773`*^9}, {3.6902687324695253`*^9, 
  3.6902687815005054`*^9}, {3.6902688306814895`*^9, 3.6902689058589983`*^9}, {
  3.6902689387096553`*^9, 3.6902690344515934`*^9}, {3.6902740259788923`*^9, 
  3.6902740586795464`*^9}, {3.6902741308830156`*^9, 3.6902741353831053`*^9}, {
  3.6902741929262867`*^9, 3.6902742030484905`*^9}, {3.6902743226689625`*^9, 
  3.6902743783520803`*^9}}],

Cell[BoxData[
 DynamicModuleBox[{$CellContext`held$$ = Hold[
   ListLinePlot[{$CellContext`XYDB, $CellContext`XYMADX}, 
    PlotStyle -> {{Blue}, {Red, Dashed}}, 
    PlotLegends -> {"Survey Database", "Survey MADX"}, 
    PlotRange -> {{1831, 1882}, {2075, 2126}}, ImageSize -> 1000, AspectRatio -> 
    1]], $CellContext`head$$ = ListLinePlot, $CellContext`plot$$ = Legended[
   Graphics[{{}, {{{}, {}, {
        Hue[0.67, 0.6, 0.6], 
        Directive[
         RGBColor[0, 0, 1]], 
        Line[CompressedData["
1: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
BYfXDw9z9D9SqCs1
         "]]}, {
        Hue[0.9060679774997897, 0.6, 0.6], 
        Directive[
         RGBColor[1, 0, 0], 
         Dashing[{Small, Small}]], 
        Line[CompressedData["
1:eJyt2Xk4Vfv3B3DFTYOUkAwZQoloIjOLzFN0KglRGUoUCZVkKCLKnGl/aqco
JEkDEhpNt2hSlG5JA+KiMobvOv78Oc+j+8t5zj+v55w/zt6ftdd6r+dIbN+7
3nkqG76mjr3ZXO4ujGbfTkDla2fD3wdo+L9ueW3E7rWNQHbO0GDgQRq81Y6F
1DkSeDG0VU7v0Hj71rlfWo7eq9tdOd2fhqwhdl5qK4GFZo8Few+Pt3SiUvqI
PYEbVt8ibwbQoFJjtn9gC4FX8xWiYoNoOHqlyZtpyRkhLL2vP9BLlvn9uQ+7
5YNpaJbJr59tQyC87mwvRwgNp+45B2zcTEA1f5lmzCQ4wrbakbImEJVsJCBy
lIZLmwRzczYSWFHrbRN/jIYPEm1zw9AOUVH389A/injmOmwg0Pku7vypUBpc
A0IzVdFn+d2+ZKNf5d6ctZBBwPKl3rfgMBrYl5hcLFtPQNtVwL+ShWN26oeZ
o1WEoq+LHh9vI938UKazOUfHDLL8g6oWBEZrfmnSERO7+NALt5vmBArqv6U2
szAn24s1F0wJGD7L4e48Md4+H4/H2qP9/APzeSLHm5fPro/HhMAzrs7eQvQF
owM5KUYEHqzuqfCNoqFi9NijcgMCHd/EGM4nadCrOtH3SIdAab/fnf3RNHQ2
f+xlmlfvKkv72BqzxwCBDPm1Z5vQIqEcwS1aBAKrpJ96xNDgtHv16VpNAqfV
BpVrJsGeHcvCizQIHAzhnLI0loa0uw89alUJXInaeUkqjoZn3+Lsi1QIzDm0
cl/IJDjNhEo7r0xgCW+SxHt05D5rwX4lAmx5Bxme8TTkn3vtWK5I4Hiwj18D
C2sLDjowvSWmfMwxOwx1eJYTiOV4eGR74sQ2SdDau0ieQHNyl041epp8/icu
OQLncsIk1U7TcHLmuqKcpfi5dsmqY/8P77FP+it8CYHpRbe4+li4p/DULENp
/P0uGjM9k2ionHMsXEiMwD0xm/gdyfh8PSo9zrToag+WfnP7xUD2QgJCVjcG
wtFUYX5GvjCBqf88dSpBh2oZWrkJEVB2tyhomgTHXltzUUIQ68H7lNQomqvI
vySPn0CCsv6HzSk0vPthONuVj4DRvFmMQ7/hzBKxO6K8BJ5KSN9OQ7sW/dop
zUOgKE9ncS2aAuue63MIzJdRufCDhXMVLbuZ3t2+fsy6el8eGc7A8825wDBO
ndge+2I//81JwCYvf6nJJNh2G5tmFzuBQgPYbPYbLvqpKZw9QkHYhmmMdeje
W4Eh6v0UvN7k8WUKum9+cTDTzxdtYWmRiH0WNb0U6B0qFX+H169z7ef8pz8o
8Fsg2n0ZvUVrWD7iOwVrHK/ynZwExzyy7tbuoeBHjI6AB9ODLw7+3UlBx4c2
bQG0gLKf0tEOClZwlej0JU9snuPN69W+UVAVNtv8FVrpqrMwtFHwS8pjYSY6
41Dml5RWCnI0rPoSWZiXbXDkyVcKtmyShRAWTgnoGWZavWTumKs+DWmeaabg
fbdiQWLSxLbTbeLMe0/BWbM0txUs/LagRPzfdxQIvuWofYHPs+i8ymvLGimY
zXct7Q32k2HC1lj9moLNbnc14tFOxYcb+J5TcH3z44y5CTTUvHzxmmmb6kSW
3nTFgvvtUwpsp7hKZWF/s+JO4VCoo+DwtdF6IXTvyeM+ArUUxDsVZJXG/bm5
NxreH31MwW03LpHtaL68Li2eGgq0eBXv7cJ50HRPneNoNQWlV7Py1X7Dun3s
ks2VFLgYtJVexfmiYW0hqFiB59cb4rUIvfVGUfjgQwrmuVRfKY4e71VvG44z
Xc0IGvPw2yafvnKsl5hrX6+enNido2sc2srw/r/nj5KbBKfvbRWchW6tTG0p
wfne9lT5peQdCvJN3VxqcP4/Sy+MTi6hwOeAtWMYWugvTtP8Ygo0hhY3J2Ce
6Ck/bepRREFQa5vFD8wfJgKffS1uUOA4c4+rH+abs0E+Y1ZKEmJp0ch1BYev
UxDoYN7xBvNTqsSJdXYFFGRwyC2bhq7YEh107hoFSTKVN/xD/9y8no/+/ZRP
gbFFWV4X5r/DnBznOvOwXl5GpvtjXuQ3EudajZYQPpDSFfLn7i/OvOJ3hYKr
Cn5Fzmg5xgcT61wKhp5t6I/GPDzzROvU0MsUJPjyi/wIGu/Lbowx267KFGba
qsn4x74sChQ4zKNEjkzsGuU9fUbogZaTAW2Y3z3NfOqDMygIkLrlvwv3hdeD
Nk3W6Cntnl7KaLuXYjWe57Ee8tuXS/hi/rvSV800H3cES982io03RP/Ls7H9
ts94exm7jkijtbgSdu5Cq3yXlylLp+BcTEip+34azHRLThxG57hJ81Z5/7lH
A69fUkVb8+uXSaFfxRksDztHQUv6l6rdXtiP9XVN16KbX0l7VXhO7Ef7wvzZ
0Lsvq5pIoqNEzDNSaApqV8jQi/fSYMHhmGWBvnHgpdOePePd4S875tIeqzHP
1w0zyTiL5390Ya3y7onNbnkx4zxa6r5kXIDbf3epZ7bDOXRZU8yMw7sm9pNZ
UVVn0d0/7jAqd9LgLiThTNDZAlM3FLnifPzk0Euh+/3DrMLR54mrexr6XG73
UIPLeEfEHXVlmv3Q16evWFi56K5pCjq+yWjgqzMNlgM9ixTxfmndaL2cvx3n
9TZBSabLWsJY+uuw/iVrtGXhXKNN6F2i4Qcc0E5CvvK66NDjllf90duFPsYN
bMO8Oi92zjt09yUpv1+ONJyI/6dyIZ7vnMHA4M2T4KaXxvft0aqXN7YUONBg
WNfrNQPrsU09QWbQHvNSVMUeY7RTi0boxt9wv1q+dQTadEiIP9+OBpkH3Tfb
0Y0RdgMFtjSsVPHzWo3PV0Du4t1DW8Z7r3v2mBd2pLsxXW62LZILn/fC0Isz
P1pP7N3G67ZFoOUrNL1NfsOewf1lLZkUrHpQ+vrcRhpom9PRjhcpCMnrC56J
VronZ69yiYIjMUuKpm74c69SNzzAg/0tKjNr88B63FeKFDLnY7/sDnfdLGqJ
+9cLMuagJF6W9tYduvQJXeZX1xO8joba/Kq1U7FfuyX97Mi2oCGrzeBcMdpD
TvOC2CSYk+1qkTfOB7lMs8fx5jQI2zh6OOE8ql8jESlihnkhQf7UQpxXax40
5MaYTmy1gu/e9Wjvh7LDHOh34ve4v+P8FP20rl/ehIbrYqsT7uF8lQn/a0Ws
8Xjr6VWMOT4tcTnTWirbxesLKShwTMqcYzixV2Z1MdZhPkl09o41WUvDl67W
RcaYb1IrhZao69Cg+SpGkmm9gObFrEwpyAkEPaLgzN51+s+Ahpz7e1PkMF9l
n3fn1UU3C8uGQxXOE3ODnBbtP7fXX4ZxmzDPieWc+BiGXlaaexkwH743XaoR
qYX18LSqfNMTCm4e4vaS/w3r/1N/xB3zZn3o4L1aTRpm+wpPWYD5dlqRQKw2
+u/5qUXrnlFgeF/Yo15jvHvpA8VM779e5860brpsw8grClJSapzN1Sf2z5wA
+/kNFBxSyKic9xv+5btom1MTBfPvPs+XU/vv5t3azF31gQJaiNGorUpDXNK3
2J9fKHhxWcnKTIWG7+/z45jekCrB0grRdVf4cX+51nCTy5iFLRvTxJ1xH9JZ
39C2GB3ysHT/1nYKuiJKAr8r07D2Y+XjDNyvVF0qVjybBNtcUfXuxP3NKWe0
7Co6/cgPe8duzG8hRuzW6EOfK3iycD/8NndmrNIkWMTO0r4b98+l5juX8qEt
DL3z237i/hRqWvp+DQ3tVxbI2PZRINyyAW6zcGOx2FKm28TUxnwwNJxROYz9
9cKXGJffcP3O54GDuJ/3JXV8dEZLCz6Z8hn3d7UhG5Ed6F2dhdEy0wisXx3w
ctsk2OTZylHn6QQsHx2LdWR+vqor6QU3gX84Fw7nonc+v5HMdIv6CpYu8hyu
0ZtLwK645v5TdLrDIXbteQTSmnd/n4r3T/ciX3cYLwHGZbO6RZPgWpcFS+v4
CBwXNXqsi2ZsSaPXLiDAvevMhWT0zx+prVGCBHQs4XzhJJiz/WTFSyECvgf9
k16jVV7PHK0VIXDznHDFDKz/8BRfYR1R/L7B1gvyLGwidESEacf3fmOeLe7z
1kuKgHXv+pVCqhPbhlf6uLY0Ac5bd3X00Bq6RxP95Anotxon9eDzrszwAn0F
AiJbH14Uwf7B9fknDK0m0FdVqHUW+1VHJdFh2te6i6UDxKW4FJUIaAyLlk/D
flgzpZG/eQ0BaeGw2Ylov1TdnH5lAssY7Rrdk2CNOcus56oS2CB5sNIc+/MH
qaNvRtQJtEtZLTLAfl9XPSN5viYBU/utBem/4dRVnn3yWgRWXUuVHkV/5V08
5a02Xn/dFO+TOJ9kWopj+XUIHNYsVvrJwlxlqXFM05Yca5jelJN0YYEBgTtB
T4dGdCf2fMfQwPvog5m8WYy1E9su2TQk0ZCAQrJS6W1030lHMDInUPVY+fop
A8x3/u+fpKPFzJITb6PLDOPO37UiMHyNL0AF57/bu6kXmJZZvJOldVea2g+h
q6MYBjtY2FZoqWv8egJP9ZX6nqP5WjqrCxhYX0HWhVmYR+TKOM7s2UCgbPHD
lYqmf+4SvgJZ2Y0EPGfdVyxD71zie1DCmoCPaa+OIuYrOzfJL87oZ0VzBl1Y
+LvCner96DBxFfZ96NwV5V3BmwlsV90WrIp5jfv7893JaKeKmXIm6IqSlNpC
GwIpU+S/iWFenPskViByC4EefffK2yx8q33JAqanR74Zs+2jVYuJPYGpDoqb
49dPbK+3m1a0ok3Us0Y7WNiya+vorK143qdtAwUYNGTca58y6kDgvr+Jmxvm
7QUBNfv3OBIwup5n9O9vmHN9490RdFhgjlXaJpzPIstvRG8jsHdts7gR5nvF
j1KahehqtrneB9B+Xo7rzu4gwDU18fop3C8UksPH/CGjj6VjrARF2tDiZ7/5
G+D+0iLE6LF0wn7HPSK9FvebyMEbMfHooGk5zz9PghMkH3vVo1381aUjcZ+6
06DeHuGM/aghWUgP97O0jVbPKXRTza5T+9CyDJ7Hv9CG4u12W3Gf04/+dWaO
CwHNcoPUKHTzpjU+vui8FZv1InCfnPZQ0usuetX77ECd7eO909xmzEtdFIKY
Pjw6UuvmSkD27nb+A84TW0WMbb87+mEHd9dKl/9u81brAQ90uaWt6jv0Kx6r
jZ5oJyWDw63MfZ57/0tv9OqVbzl9d03sFz0Pn/iib7kotz9yo0Fn6em+N3h9
qvRp5REPGjrv7R1zkecAS6v9q+KWgM5dS6xu7sH80dA/cxXauNzWo2AvDTtW
P7Dswfv/0cbK2drzz70sB6yvoVNrZNl+oe8fTLg0DV1r8e32jX00zGiMyhzA
+uA3Gl79Cn3YViGPWT/TG+Wim7zx+jIb74Uxfca5dgS9/1TqDin0BSGp9qH9
NFh58me9xXpOHrU5eN1nvO3Umi8xPVDWNOaIB8vG/r/3t28Y+/9+Iv8PRtM4
FA==
         "]]}}}, {}}, {
    AxesLabel -> {None, None}, AspectRatio -> 1, Axes -> True, 
     AxesOrigin -> {1840., 2076.}, ImageSize -> 1000, Method -> {}, 
     PlotRange -> {{1831, 1882}, {2075, 2126}}, PlotRangeClipping -> True, 
     PlotRangePadding -> {Automatic, Automatic}}], {
    Placed[
     LineLegend[{
       Directive[
        RGBColor[0, 0, 1]], 
       Directive[
        RGBColor[1, 0, 0], 
        Dashing[{Small, Small}]]}, {"Survey Database", "Survey MADX"}, 
      LegendLayout -> "Column", LegendMarkers -> False], After, 
     Identity]}], $CellContext`range$$ = {{1831., 1882.}, {2075., 
  2126.}}, $CellContext`zoomRange$$ = {{1831., 1882.}, {2075., 
  2126.}}, $CellContext`defRange$$ = {{1831., 1882.}, {2075., 
  2126.}}, $CellContext`defSize$$ = {1137, 
  1002}, $CellContext`defEpilog$$ = {}, $CellContext`size$$ = {1137, 
  1002}, $CellContext`pt1$$ = {1855.6705202241096`, 
  2125.1677534020105`}, $CellContext`pt2$$ = {1855.6695414897883`, 
  2125.169539498398}, $CellContext`spt1$$ = {0.34800000000000003`, 
  0.8560000000298025}, $CellContext`dist$$, $CellContext`rect$$ = {}, \
$CellContext`dragged$$ = False, $CellContext`overButton$$ = 
  False, $CellContext`pushedButton$$ = False, $CellContext`overHor$$ = 
  False, $CellContext`overVer$$ = False, $CellContext`activeHor$$ = 
  False, $CellContext`activeVer$$ = False, $CellContext`xPos$$ = 
  0.5, $CellContext`yPos$$ = 0.5, $CellContext`minDistance$$ = 
  0.1, $CellContext`sliderRatio$$ = 0.05, $CellContext`panFactor$$ = 
  10, $CellContext`reset$$, $CellContext`distance$$, $CellContext`slider$$, \
$CellContext`replot$$, $CellContext`coord$$ = None}, 
  TagBox[
   TagBox[
    DynamicBox[ToBoxes[
      MouseAppearance[
       Show[$CellContext`plot$$, PlotRange -> Dynamic[$CellContext`range$$], 
        Epilog -> {$CellContext`defEpilog$$, 
          EventHandler[
           Inset[
            Graphics[{Black, 
              Dynamic[
               Opacity[
                If[$CellContext`overButton$$, 
                 If[$CellContext`pushedButton$$, 0.2, 0.1], 0]]], 
              Rectangle[{0, 0}, {1, 1}, RoundingRadius -> Offset[5]], Black, 
              Dynamic[
               Opacity[
                If[$CellContext`overButton$$, 
                 If[$CellContext`pushedButton$$, 0.9, 0.7], 0]]], 
              Text["Replot", Center]}, ImageSize -> 50, AspectRatio -> 1/2, 
             BaseStyle -> {FontFamily -> "Helvetica"}], 
            Scaled[{1, 1}], 
            Scaled[{1, 1}]], {
           "MouseEntered" :> ($CellContext`overButton$$ = True), 
            "MouseExited" :> ($CellContext`overButton$$ = False), 
            "MouseDown" :> ($CellContext`pushedButton$$ = True), 
            "MouseUp" :> ($CellContext`pushedButton$$ = 
              False; $CellContext`replot$$[$CellContext`plot$$])}], 
          FaceForm[{Blue, 
            Opacity[0.1]}], 
          EdgeForm[{Thin, Dotted, 
            Opacity[0.5]}], 
          Dynamic[$CellContext`rect$$], 
          Inset[
           $CellContext`slider$$[
            Dynamic[$CellContext`xPos$$], 
            Dynamic[$CellContext`overHor$$], 
            Dynamic[$CellContext`activeHor$$], $CellContext`size$$, True], 
           Scaled[{0.5, 0}], 
           Scaled[{0.5, 0}]], 
          Inset[
           $CellContext`slider$$[
            Dynamic[$CellContext`yPos$$], 
            Dynamic[$CellContext`overVer$$], 
            Dynamic[$CellContext`activeVer$$], $CellContext`size$$, False], 
           Scaled[{0, 0.5}], 
           Scaled[{0, 0.5}]]}], 
       Which[$CellContext`dragged$$, "FrameRisingResize", 
        CurrentValue["AltKey"], 
        Graphics[{Antialiasing -> False, 
          Line[{{-1, 0}, {1, 0}}], 
          Line[{{0, -1}, {0, 1}}], 
          Dynamic[
           Text[
            Style[$CellContext`coord$$, FontFamily -> "Helvetica", 8], {0, 
            0}, {-1.1, -1}]]}, ImageSize -> 160, ImagePadding -> 0, 
         AspectRatio -> 1, PlotRange -> {{-1, 1}, {-1, 1}} 8], 
        And[
         CurrentValue["ControlKey"], 
         Not[$CellContext`overHor$$], 
         Not[$CellContext`overVer$$]], "ZoomView", 
        And[
         CurrentValue["ShiftKey"], 
         Not[$CellContext`overHor$$], 
         Not[$CellContext`overVer$$]], "DragGraphics", 
        And[
         CurrentValue["ShiftKey"], 
         Or[$CellContext`overHor$$, $CellContext`overVer$$]], "LinkHand", 
        True, Automatic]], StandardForm],
     ImageSizeCache->{1125., {498., 502.}},
     TrackedSymbols:>{$CellContext`dragged$$, $CellContext`plot$$}],
    EventHandlerTag[{
     "MouseMoved" :> ($CellContext`coord$$ = MousePosition["Graphics"]), 
      "MouseClicked" :> If[
        And[
         Not[$CellContext`overButton$$], CurrentValue["MouseClickCount"] == 
         2], 
        $CellContext`reset$$[]], "MouseDown" :> If[
        Or[
         Not[$CellContext`overHor$$], 
         Not[$CellContext`overVer$$]], $CellContext`pt1$$ = 
         MousePosition["Graphics"]; $CellContext`spt1$$ = 
         MousePosition["GraphicsScaled"]], "MouseDragged" :> If[
        And[
         Or[
          And[
           Not[$CellContext`overHor$$], 
           Not[$CellContext`overVer$$]], $CellContext`dragged$$], 
         Not[$CellContext`overButton$$], 
         Not[$CellContext`pushedButton$$], 
         Not[$CellContext`activeVer$$], 
         Not[$CellContext`activeHor$$]], $CellContext`dragged$$ = True; Which[
          CurrentValue["ControlKey"], $CellContext`dist$$ = Last[
              MousePosition["GraphicsScaled"]] - 
            Last[$CellContext`spt1$$]; $CellContext`range$$ = \
$CellContext`zoomRange$$ 10^$CellContext`dist$$, 
          CurrentValue["ShiftKey"], $CellContext`pt2$$ = MapThread[Rescale, {
               MousePosition["GraphicsScaled"], {{0, 1}, {0, 
               1}}, $CellContext`zoomRange$$}] - $CellContext`pt1$$; \
$CellContext`range$$ = $CellContext`zoomRange$$ - $CellContext`pt2$$, 
          True, $CellContext`pt2$$ = 
           MousePosition["Graphics"]; $CellContext`rect$$ = 
           If[$CellContext`pt2$$ === None, {}, 
             Rectangle[$CellContext`pt1$$, $CellContext`pt2$$]]]], 
      "MouseUp" :> (If[
         And[
          Not[
           CurrentValue["ShiftKey"]], 
          Not[
           CurrentValue[
           "ControlKey"]], $CellContext`dragged$$, \
$CellContext`distance$$[$CellContext`pt1$$, $CellContext`pt2$$, \
$CellContext`zoomRange$$] > $CellContext`minDistance$$], $CellContext`rect$$ = \
{}; $CellContext`range$$ = Transpose[
            Sort[{$CellContext`pt1$$, $CellContext`pt2$$}]]; 
         Null]; $CellContext`zoomRange$$ = $CellContext`range$$; \
$CellContext`dragged$$ = False), PassEventsDown -> True, PassEventsUp -> True,
       PassEventsDown -> Automatic, PassEventsUp -> True}]],
   Deploy,
   DefaultBaseStyle->"Deploy"],
  DynamicModuleValues:>{{DownValues[$CellContext`reset$$] = {HoldPattern[
         $CellContext`reset$$[]] :> ($CellContext`range$$ = \
($CellContext`zoomRange$$ = $CellContext`defRange$$); $CellContext`xPos$$ = \
($CellContext`yPos$$ = 0.5))}}, {
    DownValues[$CellContext`distance$$] = {HoldPattern[
         $CellContext`distance$$[
          Pattern[$CellContext`p1, 
           Blank[]], 
          Pattern[$CellContext`p2, 
           Blank[]], {
           Pattern[$CellContext`r1, 
            Blank[]], 
           Pattern[$CellContext`r2, 
            Blank[]]}]] :> Quiet[
         N[
          EuclideanDistance[
           Rescale[$CellContext`p1, $CellContext`r1], 
           Rescale[$CellContext`p2, $CellContext`r2]]]]}}, {
    DownValues[$CellContext`slider$$] = {HoldPattern[
         $CellContext`slider$$[
          Dynamic[
           Pattern[$CellContext`pos$, 
            Blank[]]], 
          Dynamic[
           Pattern[$CellContext`over$, 
            Blank[]]], 
          Dynamic[
           Pattern[$CellContext`active$, 
            Blank[]]], 
          Pattern[$CellContext`size$, 
           Blank[]], 
          Pattern[$CellContext`hor$, 
           Blank[]]]] :> EventHandler[
         Dynamic[
          LocatorPane[
           If[$CellContext`hor$, 
            Dynamic[{$CellContext`pos$, 0.5}, ($CellContext`pos$ = Clip[
                First[#], {0, 1}]; Part[$CellContext`range$$, 1] = If[
                CurrentValue["ShiftKey"], 
                First[$CellContext`zoomRange$$] + ($CellContext`pos$ - 0.5) 
                 Abs[
                   Apply[Subtract, 
                    
                    First[$CellContext`zoomRange$$]]] \
$CellContext`panFactor$$, First[$CellContext`zoomRange$$] 
                10^(($CellContext`pos$ - 0.5) 2)])& ], 
            Dynamic[{0.5, $CellContext`pos$}, ($CellContext`pos$ = Clip[
                Last[#], {0, 1}]; Part[$CellContext`range$$, 2] = If[
                CurrentValue["ShiftKey"], 
                Last[$CellContext`zoomRange$$] + ($CellContext`pos$ - 0.5) Abs[
                   Apply[Subtract, 
                    Last[$CellContext`zoomRange$$]]] $CellContext`panFactor$$,
                 Last[$CellContext`zoomRange$$] 
                10^(($CellContext`pos$ - 0.5) 2)])& ]], 
           Graphics[{Black, 
             Opacity[
              If[
               And[
                Or[$CellContext`over$, $CellContext`active$], 
                Not[$CellContext`dragged$$]], 0.1, 0]], 
             Rectangle[{0, 0}, {1, 1}, RoundingRadius -> Offset[5]], 
             Opacity[1], 
             Locator[
              Dynamic[
               
               If[$CellContext`hor$, {$CellContext`pos$, 0.5}, {
                0.5, $CellContext`pos$}]], Appearance -> If[
                And[
                 Or[$CellContext`over$, $CellContext`active$], 
                 Not[$CellContext`dragged$$]], Automatic, None]]}, 
            ImagePadding -> 0, PlotRangePadding -> 0, ImageSize -> 
            If[$CellContext`hor$, # {1, 0} + {0, 20}& , # {0, 1} + {20, 
               0}& ][$CellContext`size$], AspectRatio -> 
            If[$CellContext`hor$, (1/$CellContext`sliderRatio$$)/
              First[$CellContext`size$], $CellContext`sliderRatio$$ 
              Last[$CellContext`size$]]], Enabled -> 
           Or[$CellContext`over$, $CellContext`active$], Appearance -> None, 
           FrameMargins -> 0, ContentPadding -> False], 
          TrackedSymbols :> {$CellContext`pos$, $CellContext`over$, \
$CellContext`active$, $CellContext`dragged$$}], {
         "MouseDown" :> ($CellContext`active$ = True), 
          "MouseUp" :> ($CellContext`active$ = False), 
          "MouseEntered" :> ($CellContext`over$ = Not[
             
             If[$CellContext`hor$, $CellContext`activeVer$$, \
$CellContext`activeHor$$]]), "MouseExited" :> ($CellContext`over$ = False)}, 
         PassEventsDown -> True]}}, {
    DownValues[$CellContext`replot$$] = {HoldPattern[
         $CellContext`replot$$[
          Pattern[$CellContext`p$, 
           Blank[]]]] :> Module[{$CellContext`temp$}, 
         Switch[$CellContext`head$$, 
          Plot, $CellContext`temp$ = ReplacePart[$CellContext`held$$, {{1, 2} -> {
                Part[$CellContext`held$$, 1, 2, 1], 
                Part[$CellContext`range$$, 1, 1], 
                Part[$CellContext`range$$, 1, 2]}}]; $CellContext`plot$$ = 
           ReleaseHold[
             If[
              MemberQ[$CellContext`temp$, PlotRange, Infinity], 
              ReplaceAll[$CellContext`temp$, {Blank[][PlotRange, 
                  Blank[]] -> (PlotRange -> $CellContext`range$$)}], 
              
              Insert[$CellContext`temp$, PlotRange -> $CellContext`range$$, {
               1, -1}]]], 
          Alternatives[
          DensityPlot, ContourPlot, RegionPlot, StreamPlot, StreamDensityPlot,
            VectorPlot, VectorDensityPlot, 
           LineIntegralConvolutionPlot], $CellContext`temp$ = 
           ReplacePart[$CellContext`held$$, {{1, 2} -> {
                Part[$CellContext`held$$, 1, 2, 1], 
                Part[$CellContext`range$$, 1, 1], 
                Part[$CellContext`range$$, 1, 2]}, {1, 3} -> {
                Part[$CellContext`held$$, 1, 3, 1], 
                Part[$CellContext`range$$, 2, 1], 
                Part[$CellContext`range$$, 2, 2]}}]; $CellContext`plot$$ = 
           ReleaseHold[
             If[
              MemberQ[$CellContext`temp$, PlotRange, Infinity], 
              ReplaceAll[$CellContext`temp$, {Blank[][PlotRange, 
                  Blank[]] -> (PlotRange -> $CellContext`range$$)}], 
              
              Insert[$CellContext`temp$, PlotRange -> $CellContext`range$$, {
               1, -1}]]]; Null, 
          Blank[], Null]]}}},
  Initialization:>($CellContext`slider$$[
      Dynamic[
       Pattern[$CellContext`pos$, 
        Blank[]]], 
      Dynamic[
       Pattern[$CellContext`over$, 
        Blank[]]], 
      Dynamic[
       Pattern[$CellContext`active$, 
        Blank[]]], 
      Pattern[$CellContext`size$, 
       Blank[]], 
      Pattern[$CellContext`hor$, 
       Blank[]]] := EventHandler[
      Dynamic[
       LocatorPane[
        If[$CellContext`hor$, 
         Dynamic[{$CellContext`pos$, 0.5}, ($CellContext`pos$ = Clip[
             First[#], {0, 1}]; Part[$CellContext`range$$, 1] = If[
             CurrentValue["ShiftKey"], 
             First[$CellContext`zoomRange$$] + ($CellContext`pos$ - 0.5) Abs[
                Apply[Subtract, 
                 First[$CellContext`zoomRange$$]]] $CellContext`panFactor$$, 
             First[$CellContext`zoomRange$$] 
             10^(($CellContext`pos$ - 0.5) 2)])& ], 
         Dynamic[{0.5, $CellContext`pos$}, ($CellContext`pos$ = Clip[
             Last[#], {0, 1}]; Part[$CellContext`range$$, 2] = If[
             CurrentValue["ShiftKey"], 
             Last[$CellContext`zoomRange$$] + ($CellContext`pos$ - 0.5) Abs[
                Apply[Subtract, 
                 Last[$CellContext`zoomRange$$]]] $CellContext`panFactor$$, 
             Last[$CellContext`zoomRange$$] 
             10^(($CellContext`pos$ - 0.5) 2)])& ]], 
        Graphics[{Black, 
          Opacity[
           If[
            And[
             Or[$CellContext`over$, $CellContext`active$], 
             Not[$CellContext`dragged$$]], 0.1, 0]], 
          Rectangle[{0, 0}, {1, 1}, RoundingRadius -> Offset[5]], 
          Opacity[1], 
          Locator[
           Dynamic[
            
            If[$CellContext`hor$, {$CellContext`pos$, 0.5}, {
             0.5, $CellContext`pos$}]], Appearance -> If[
             And[
              Or[$CellContext`over$, $CellContext`active$], 
              Not[$CellContext`dragged$$]], Automatic, None]]}, ImagePadding -> 
         0, PlotRangePadding -> 0, ImageSize -> 
         If[$CellContext`hor$, # {1, 0} + {0, 20}& , # {0, 1} + {20, 
            0}& ][$CellContext`size$], AspectRatio -> 
         If[$CellContext`hor$, (1/$CellContext`sliderRatio$$)/
           First[$CellContext`size$], $CellContext`sliderRatio$$ 
           Last[$CellContext`size$]]], Enabled -> 
        Or[$CellContext`over$, $CellContext`active$], Appearance -> None, 
        FrameMargins -> 0, ContentPadding -> False], 
       TrackedSymbols :> {$CellContext`pos$, $CellContext`over$, \
$CellContext`active$, $CellContext`dragged$$}], {
      "MouseDown" :> ($CellContext`active$ = True), 
       "MouseUp" :> ($CellContext`active$ = False), 
       "MouseEntered" :> ($CellContext`over$ = Not[
          
          If[$CellContext`hor$, $CellContext`activeVer$$, \
$CellContext`activeHor$$]]), "MouseExited" :> ($CellContext`over$ = False)}, 
      PassEventsDown -> True]; $CellContext`replot$$[
      Pattern[$CellContext`p$, 
       Blank[]]] := Module[{$CellContext`temp$}, 
      Switch[$CellContext`head$$, 
       Plot, $CellContext`temp$ = ReplacePart[$CellContext`held$$, {{1, 2} -> {
             Part[$CellContext`held$$, 1, 2, 1], 
             Part[$CellContext`range$$, 1, 1], 
             Part[$CellContext`range$$, 1, 2]}}]; $CellContext`plot$$ = 
        ReleaseHold[
          If[
           MemberQ[$CellContext`temp$, PlotRange, Infinity], 
           ReplaceAll[$CellContext`temp$, {Blank[][PlotRange, 
               Blank[]] -> (PlotRange -> $CellContext`range$$)}], 
           
           Insert[$CellContext`temp$, PlotRange -> $CellContext`range$$, {
            1, -1}]]], 
       Alternatives[
       DensityPlot, ContourPlot, RegionPlot, StreamPlot, StreamDensityPlot, 
        VectorPlot, VectorDensityPlot, 
        LineIntegralConvolutionPlot], $CellContext`temp$ = 
        ReplacePart[$CellContext`held$$, {{1, 2} -> {
             Part[$CellContext`held$$, 1, 2, 1], 
             Part[$CellContext`range$$, 1, 1], 
             Part[$CellContext`range$$, 1, 2]}, {1, 3} -> {
             Part[$CellContext`held$$, 1, 3, 1], 
             Part[$CellContext`range$$, 2, 1], 
             Part[$CellContext`range$$, 2, 2]}}]; $CellContext`plot$$ = 
        ReleaseHold[
          If[
           MemberQ[$CellContext`temp$, PlotRange, Infinity], 
           ReplaceAll[$CellContext`temp$, {Blank[][PlotRange, 
               Blank[]] -> (PlotRange -> $CellContext`range$$)}], 
           
           Insert[$CellContext`temp$, PlotRange -> $CellContext`range$$, {
            1, -1}]]]; Null, 
       Blank[], 
       Null]]; $CellContext`reset$$[] := ($CellContext`range$$ = \
($CellContext`zoomRange$$ = $CellContext`defRange$$); $CellContext`xPos$$ = \
($CellContext`yPos$$ = 0.5)); $CellContext`distance$$[
      Pattern[$CellContext`p1, 
       Blank[]], 
      Pattern[$CellContext`p2, 
       Blank[]], {
       Pattern[$CellContext`r1, 
        Blank[]], 
       Pattern[$CellContext`r2, 
        Blank[]]}] := Quiet[
      N[
       EuclideanDistance[
        Rescale[$CellContext`p1, $CellContext`r1], 
        Rescale[$CellContext`p2, $CellContext`r2]]]]; Null)]], "Output",
 CellChangeTimes->{
  3.690265600146309*^9, 3.690265725759925*^9, 3.690265757303731*^9, 
   3.6902658057894635`*^9, {3.6902659207012734`*^9, 3.690265940841261*^9}, 
   3.6902660323058195`*^9, {3.6902660914607954`*^9, 3.6902660960628843`*^9}, 
   3.690267454684072*^9, 3.690267488670761*^9, {3.69026808178866*^9, 
   3.6902681385288315`*^9}, 3.690268185512792*^9, 3.6902682318277216`*^9, 
   3.6902682914349337`*^9, 3.690268322713578*^9, 3.690268361544355*^9, 
   3.6902683945930405`*^9, {3.690268748449845*^9, 3.690268807461025*^9}, {
   3.690268838395647*^9, 3.6902689066290135`*^9}, {3.6902689415297117`*^9, 
   3.690269035741619*^9}, 3.6902740597905684`*^9, 3.6902741383281736`*^9, 
   3.6902742039485083`*^9, {3.6902743164088373`*^9, 3.690274326079031*^9}, 
   3.690274383292179*^9, 3.6902744756322236`*^9}]
}, Open  ]],

Cell[BoxData[
 RowBox[{"\[AliasDelimiter]", "\[AliasDelimiter]"}]], "Input",
 CellChangeTimes->{{3.690268799030856*^9, 3.690268799910874*^9}}],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"2", " ", "\[Pi]", " ", 
  RowBox[{"200", "/", "16"}]}]], "Input",
 CellChangeTimes->{{3.690274630661561*^9, 3.690274670504423*^9}, {
  3.690274836652116*^9, 3.690274857682536*^9}}],

Cell[BoxData[
 RowBox[{"25", " ", "\[Pi]"}]], "Output",
 CellChangeTimes->{3.690274858584556*^9}]
}, Open  ]]
}, Open  ]]
},
WindowSize->{1305, 1012},
WindowMargins->{{-1, Automatic}, {Automatic, 0}},
PrivateNotebookOptions->{"VersionedStylesheet"->{"Default.nb"[8.] -> True}},
FrontEndVersion->"9.0 for Microsoft Windows (32-bit) (January 25, 2013)",
StyleDefinitions->"Default.nb"
]
(* End of Notebook Content *)

(* Internal cache information *)
(*CellTagsOutline
CellTagsIndex->{}
*)
(*CellTagsIndex
CellTagsIndex->{}
*)
(*NotebookFileOutline
Notebook[{
Cell[CellGroupData[{
Cell[1485, 35, 78, 1, 99, "Title"],
Cell[CellGroupData[{
Cell[1588, 40, 305, 6, 52, "Subsection"],
Cell[1896, 48, 1015, 23, 132, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[2948, 76, 347, 6, 52, "Subsection"],
Cell[3298, 84, 1119, 25, 132, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[4454, 114, 364, 5, 52, "Subsection"],
Cell[4821, 121, 33730, 747, 1252, "Input"]
}, Open  ]]
}, Closed]],
Cell[CellGroupData[{
Cell[38600, 874, 956, 17, 181, "Title"],
Cell[39559, 893, 2588, 43, 72, "Input"],
Cell[CellGroupData[{
Cell[42172, 940, 411, 8, 31, "Input"],
Cell[42586, 950, 692, 15, 30, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[43315, 970, 259, 6, 31, "Input"],
Cell[43577, 978, 1528, 34, 50, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[45142, 1017, 232, 4, 31, "Input"],
Cell[45377, 1023, 3333, 46, 30, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[48747, 1074, 136, 2, 31, "Input"],
Cell[48886, 1078, 94, 1, 30, "Output"]
}, Open  ]],
Cell[48995, 1082, 1601, 32, 52, "Input"],
Cell[CellGroupData[{
Cell[50621, 1118, 1497, 33, 52, "Input"],
Cell[52121, 1153, 29767, 591, 1017, "Output"]
}, Open  ]],
Cell[81903, 1747, 142, 2, 31, "Input"],
Cell[CellGroupData[{
Cell[82070, 1753, 204, 4, 31, "Input"],
Cell[82277, 1759, 97, 2, 86, "Output"]
}, Open  ]]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature xwT7PERgJG8jeCg2Z9p2WINw *)
