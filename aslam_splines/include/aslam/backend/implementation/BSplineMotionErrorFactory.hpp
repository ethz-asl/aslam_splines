
namespace aslam {
namespace backend {

template<class SplineDv>
void addMotionErrorTerms(OptimizationProblemBase & problem, SplineDv & splineDv, Eigen::MatrixXd W, unsigned int errorTermOrder) {
  // Add one error term for each segment
  auto spline = splineDv.spline();

  for(int i = 0; i < spline.numValidTimeSegments(); ++i) {
    Eigen::MatrixXd R = spline.segmentIntegral(i, W, errorTermOrder);
    Eigen::VectorXd c = spline.segmentCoefficientVector(i);
    // These indices tell us the order of R and c.
    Eigen::VectorXi idxs = spline.segmentVvCoefficientVectorIndices(i);

    std::vector< DesignVariable * > dvs;
    for(unsigned i = 0; i < idxs.size(); ++i) {
      dvs.push_back(splineDv.designVariable(idxs[i]));
    }
    boost::shared_ptr< MarginalizationPriorErrorTerm > err( 
        new MarginalizationPriorErrorTerm( dvs, R*c, R ));
    problem.addErrorTerm(err);
            
  }
}


} // namespace backend
} // namespace aslam
