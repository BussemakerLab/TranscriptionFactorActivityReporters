#include<Rcpp.h>
using namespace Rcpp;
using namespace std;

class Score{
public:
  Score(Rcpp::NumericVector, double);
  Score(Rcpp::NumericVector, Rcpp::NumericVector, double);
  Rcpp::NumericMatrix scoreSeq(Rcpp::IntegerVector);
  Rcpp::NumericMatrix nuc;
  Rcpp::NumericMatrix dinuc;
  Rcpp::NumericMatrix revNuc;
  Rcpp::NumericMatrix revDinuc;
private:
  bool isDinuc;
  int k;
  double affFloor;
};

class RapidScore{
public:
  RapidScore();
  RapidScore(double);
  void add(Rcpp::NumericVector, Rcpp::Nullable<Rcpp::NumericVector>);
  void print();
  Rcpp::List scoreBulk(Rcpp::IntegerVector);
private:
  int nModels;
  double affFloor;
  std::vector<Score> models;
};

//*******************
//*** SCORE CLASS ***
//*******************
// Constructors
Score::Score(Rcpp::NumericVector nucBetas, double affinityFloor) {
  Rcpp::Function mat = Rcpp::Environment("package:base")["matrix"];
  if (nucBetas.length() % 4 != 0) {
    Rcpp::stop("Betas improperly defined!");
  }
  k = floor(nucBetas.length()/4);
  nuc = mat(nucBetas, 4, k);
  revNuc = mat(rev(nucBetas), 4, k);
  isDinuc = false;
  affFloor = affinityFloor;
}
Score::Score(Rcpp::NumericVector nucBetas, Rcpp::NumericVector dinucBetas, double affinityFloor) {
  Rcpp::Function mat = Rcpp::Environment("package:base")["matrix"];
  Rcpp::IntegerVector dinucRevIdx = Rcpp::IntegerVector::create(15,11,7,3,14,10,6,2,13,9,5,1,12,8,4,0);

  if (nucBetas.length() % 4 != 0) {
    Rcpp::stop("Betas improperly defined!");
  }
  k = floor(nucBetas.length()/4);
  if (dinucBetas.length() != 16*(k-1)) {
    Rcpp::stop("Betas improperly defined!");
  }
  nuc = mat(nucBetas, 4, k);
  revNuc = mat(rev(nucBetas), 4, k);
  dinuc = mat(dinucBetas, 16, k-1);
  // Need to manually invert the dinucleotide features
  revDinuc = mat(dinucBetas, 16, k-1);
  for (int i=0; i<k-1; i++) {   // First loop over positions
    for (int j=0; j<16; j++) {  // Next loop over dinucs
      revDinuc(j,i) = dinuc(dinucRevIdx[j], k-i-2);
    }
  }
  isDinuc = true;
  affFloor = affinityFloor;
}

// Function to score sequence
Rcpp::NumericMatrix Score::scoreSeq(Rcpp::IntegerVector seq){
  double fwd = 0;
  double rev = 0;
  Rcpp::NumericMatrix output(2, seq.length()-k+1);
  if (isDinuc) {
    for (int currWindow=0; currWindow<seq.length()-k+1; currWindow++) {
      fwd = nuc(seq[currWindow],0);
      rev = revNuc(seq[currWindow],0);
      for (int i=1; i<k; i++) {
        fwd += nuc(seq[currWindow+i], i) + dinuc(seq[currWindow+i-1]*4+seq[currWindow+i], i-1);
        rev += revNuc(seq[currWindow+i], i) + revDinuc(seq[currWindow+i-1]*4+seq[currWindow+i], i-1);
      }
      output(0, currWindow) = max(fwd, affFloor);
      output(1, currWindow) = max(rev, affFloor);
    }
  } else {
    for (int currWindow=0; currWindow<seq.length()-k+1; currWindow++) {
      fwd = nuc(seq[currWindow],0);
      rev = revNuc(seq[currWindow],0);
      for (int i=1; i<k; i++) {
        fwd += nuc(seq[currWindow+i], i);
        rev += revNuc(seq[currWindow+i], i);
      }
      output(0, currWindow) = max(fwd, affFloor);
      output(1, currWindow) = max(rev, affFloor);
    }
  }
  return(output);
}

RCPP_MODULE(scoremodule){
  Rcpp::class_<Score>( "Score" )
  .constructor<Rcpp::NumericVector, double>("constructor for nucleotide only model")
  .constructor<Rcpp::NumericVector, Rcpp::NumericVector, double>("constructor for dinucleotide model")
  .method( "scoreSeq", &Score::scoreSeq, "documentation for scoring sequence")
  // .field( "nuc", &Score::nuc, "documentation for nuc")
  // .field( "dinuc", &Score::dinuc, "documentation for nuc")
  // .field( "revNuc", &Score::revNuc, "documentation for nuc")
  // .field( "revDinuc", &Score::revDinuc, "documentation for nuc")
  ;
}

//************************
//*** RAPIDSCORE CLASS ***
//************************
// Constructor
RapidScore::RapidScore(){
  nModels = 0;
  affFloor = -1E9;
}

RapidScore::RapidScore(double affinityFloor) {
  nModels = 0;
  affFloor = affinityFloor;
}

// add model function
void RapidScore::add(Rcpp::NumericVector nucBetas, Rcpp::Nullable<Rcpp::NumericVector> dinucBetas = R_NilValue) {
  nModels++;
  // Handle dinuc betas
  if (dinucBetas.isNotNull()) {
    Rcpp::NumericVector temp(dinucBetas);
    models.push_back(Score(nucBetas, temp, affFloor));
  } else {
    models.push_back(Score(nucBetas, affFloor));
  }
}

void RapidScore::print() {
  for (int i=0; i<nModels; i++) {
    Rcpp::Rcout << "Model " << i+1 << " Nucleotide Parameters:\n";
    Rf_PrintValue(models[i].nuc);
    Rcpp::Rcout << "Dinucleotide Parameters:\n";
    Rf_PrintValue(models[i].dinuc);
  }
}

Rcpp::List RapidScore::scoreBulk(Rcpp::IntegerVector seq) {
  Rcpp::List output(nModels);

  for (int i=0; i<nModels; i++) {
    output[i] = models[i].scoreSeq(seq);
  }

  return(output);
}

RCPP_MODULE(rapidscoremodule){
  Rcpp::class_<RapidScore>( "RapidScore" )
  .constructor("documentation for default constructor")
  .constructor<double>("documentation for constructor that defines the affinity floor")
  .method( "add", &RapidScore::add, "documentation for adding a model")
  .method( "print", &RapidScore::print, "documentation for printing all betas (diagnostic)")
  .method( "scoreBulk", &RapidScore::scoreBulk, "documentation for scoring a sequence using multiple models")
  ;
}
