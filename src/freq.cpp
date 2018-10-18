#include <Rcpp.h>
#include <omp.h>
#include <regex>

using namespace Rcpp;

//' Multiply a number by two
//' 
//' @param x A single integer.
//' @param y A single integer.
//' @export
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

NumericVector ShapeCountsCpp(std::vector<int> counts, std::vector<bool> absolute, 
													bool totals = true, bool genotypic = true, bool allelic = false,
												  bool percentage = false, bool extrapolateFreq = true,
												  double minFreqAl = -1, double minFreqGt = -1) {

	if (absolute.size() == 1) {
		absolute[1] = absolute[0];
	}  

  //int const length(counts.size());
  bool returnNas = false;
  
  // calculate frequencies of alleles -----------------------------------------

  std::vector<int> countsAl(3);
  int totalAl(0);

  countsAl[0] = counts[0] * 2 + counts[1];
  countsAl[1] = counts[2] * 2 + counts[1];
  countsAl[2] = counts[3] * 2;
  totalAl = countsAl[0] + countsAl[1] + countsAl[2];

  std::vector<double> freqsAl(3);

  if (extrapolateFreq) {

  	int denumAl(countsAl[0] + countsAl[1]);
 		for (unsigned int i = 0; i < freqsAl.size() - 1; i++) {
 			freqsAl[i] = static_cast<double>(countsAl[i]) / denumAl; // extrap. freq of REF and ALT 
 		}
  	freqsAl.back() = static_cast<double>(countsAl.back()) / totalAl; // freq of missing values

  } else {

 		for (unsigned int i = 0; i < freqsAl.size(); i++) {
 			freqsAl[i] = static_cast<double>(countsAl[i]) / totalAl;
 		}

  }

  // calculate frequencies of genotypes ---------------------------------------

  int totalGt(counts[0] + counts[1] + counts[2] + counts[3]);

  std::vector<double> freqsGt(4);

  if (extrapolateFreq) {

  	int denumGt(counts[0] + counts[1] + counts[2]);
 		for (unsigned int i = 0; i < freqsGt.size() - 1; i++) {
 			freqsGt[i] = static_cast<double>(counts[i]) / denumGt; // extrap. freq of HREF, HETRO and HALT
 		}
  	freqsGt[3] = static_cast<double>(counts[3]) / totalGt; // freq of missing values
  
  } else {
 		
 		for (unsigned int i = 0; i < freqsGt.size(); i++) {
 			freqsGt[i] = static_cast<double>(counts[i]) / totalGt;
 		}

  }

  // check minimal frequencies ------------------------------------------------

	if (minFreqAl > 0) {
		if (freqsAl[1] < minFreqAl || freqsAl[2] < minFreqAl) {
			returnNas = true;
		}
	}

	if (minFreqGt > 0) {
	  if (freqsGt[1] < minFreqGt ||freqsGt[2] < minFreqGt) {
	    returnNas = true;
	  }
	}

  // build resulting frequencies ----------------------------------------------

  std::string prefixGt("");
  std::vector<double> resultGt(0);

  if (!absolute[0]) { // for genotype
  	prefixGt = "freq.";
  	resultGt = freqsGt;
 		if (percentage) {
 			prefixGt = "perc.";
 			for (unsigned int i = 0; i < resultGt.size(); i++) {
 				resultGt[i] = resultGt[i] * 100;
 			}
 		}
  } else {
  	prefixGt = "count.";
  	std::vector<double> countsDouble(counts.begin(), counts.end());
  	resultGt = countsDouble;
  	resultGt.pop_back();
  }

  std::string prefixAl("");
  std::vector< double > resultAl(0);

  if (!absolute[1]) { // for alleles
  	prefixAl = "freq.";
  	resultAl = freqsAl;
 		if (percentage) {
 			prefixAl = "perc.";
 			for (unsigned int i = 0; i < resultAl.size(); i++) {
 				resultAl[i] = resultAl[i] * 100;
 			}
 		}
  } else {
  	prefixAl = "count.";
  	std::vector<double> countsAlDouble(countsAl.begin(), countsAl.end());
  	resultAl = countsAlDouble;
  }


 	// build titles for each frequencies returned -------------------------------

 	std::vector<std::string> namesGt = {prefixGt + "gt.HOMOREF",
 																			prefixGt + "gt.HETERO",
 																			prefixGt + "gt.HOMOALT",
 																			prefixGt + "gt.MISSVAL"};

 	std::vector<std::string> namesAl = {prefixAl + "al.REF",
 																		  prefixAl + "al.ALT",
 																		  prefixAl + "al.MISSVAL"};

  // concatenate results if needed --------------------------------------------

 	if (totals) {

    resultGt.push_back(totalGt);
	 	namesGt.push_back("count.gt.TOTAL");

    resultAl.push_back(totalAl);
  	namesAl.push_back("count.al.TOTAL");

  }

  std::vector<double> resVal;
  std::vector<std::string> resNames;

  if (genotypic && allelic) {
  	
  	//resVal.reserve(resultGt.size() + resultAl.size()); // preallocate memory
    resVal.insert(resVal.end(), resultGt.begin(), resultGt.end());
    resVal.insert(resVal.end(), resultAl.begin(), resultAl.end());

    //resNames.reserve(namesGt.size() + namesGt.size()); // preallocate memory
    resNames.insert(resNames.end(), namesGt.begin(), namesGt.end());
    resNames.insert(resNames.end(), namesAl.begin(), namesAl.end());

  } else if (genotypic) {
    resVal = resultGt;
    resNames = namesGt;
  } else if (allelic) {
    resVal = resultAl;
    resNames = namesAl;
  } else {
    returnNas = true;
  }

  /*
  std::ostringstream sstream;
	sstream << returnNas;
	std::string val = sstream.str();
	Rcout << val << "\n";
  */

  NumericVector result(1, NA_REAL);
  
  if(!returnNas && resVal.size() > 0) {
  	result[0] = resVal[0];
  	if (resVal.size() > 1) {
  		for (unsigned int i = 1; i < resVal.size(); i++) {
  			result.push_back(resVal[i]);
  		}
  	}
  	result.attr("names") = resNames;
  }

  return(result);

}




//' Multiply a number by two
//' 
//' @param x A single integer.
//' @param y A single integer.
//' @export
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]




List CheapDataFrameBuilder(List a) {

    List returned_frame = clone(a);
    GenericVector sample_row = returned_frame(0);

    StringVector row_names(sample_row.length());
    for (int i = 0; i < sample_row.length(); ++i) {
        char name[5];
        //sprintf(&(name[0]), "%d", i);
        row_names(i) = name;
    }
    returned_frame.attr("row.names") = row_names;

    StringVector col_names(returned_frame.length());
    for (int j = 0; j < returned_frame.length(); ++j) {
        char name[6];
        //sprintf(&(name[0]), "X.%d", j);
        col_names(j) = name;
    }
    returned_frame.attr("names") = col_names;
    returned_frame.attr("class") = "data.frame";

    return returned_frame;

}

//' Multiply a number by two
//' 
//' @param x A single integer.
//' @param y A single integer.
//' @export
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

DataFrame CalcFreqCpp(DataFrame gt, std::vector<bool> absolute, 
                          bool totals = true, bool genotypic = true, bool allelic = false,
                          bool percentage = false, bool extrapolateFreq = true,
                          double minFreqAl = -1, double minFreqGt = -1) {

  const int nbVariants(gt.nrow());
  const int nbSamples(gt.ncol());
  const unsigned zero(0);
  const List emptyList(Dimension(1, nbSamples), zero);

  DataFrame variant(emptyList);
  CharacterVector column(nbVariants);
  CharacterVector temp(1);

  for (int i = 0; i < nbVariants; i++) {

    for (int j = 0; j < nbSamples; j++) {
      column = gt(j);
      temp(0) = column(i);
      variant(j) = temp;
    }


  }

  return(variant);

}

    /*** R
    x <- 32
    x
    */