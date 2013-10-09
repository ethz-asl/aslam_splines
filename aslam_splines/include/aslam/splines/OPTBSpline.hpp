/*
 * OPTBSpline.hpp
 *
 *  Created on: 05.08.2012
 *      Author: hannes
 */

#ifndef OPTBSPLINE_HPP_
#define OPTBSPLINE_HPP_

#include <aslam/backend/DesignVariableMappedVector.hpp>
#include <bsplines/DiffManifoldBSpline.hpp>
#include <Eigen/StdVector>
#include <boost/ptr_container/ptr_vector.hpp>
#include <aslam/backend/TransformationExpression.hpp>
#include <aslam/backend/RotationExpression.hpp>
#include <aslam/backend/EuclideanExpression.hpp>
#include <aslam/backend/FixedPointNumber.hpp>
#include <aslam/backend/GenericScalarExpression.hpp>
#include <aslam/backend/DesignVariable.hpp>
#include <aslam/backend/VectorExpression.hpp>
#include <aslam/backend/VectorExpressionNode.hpp>

#include <bsplines/manifolds/LieGroup.hpp>


namespace aslam {
	namespace splines {
		template <typename TDiffManifoldBSplineConfiguration, typename TDiffManifoldBSplineConfigurationDerived = TDiffManifoldBSplineConfiguration>
		struct DesignVariableSegmentBSplineConf : public TDiffManifoldBSplineConfigurationDerived {
			typedef TDiffManifoldBSplineConfigurationDerived ParentConf;
			typedef DesignVariableSegmentBSplineConf<TDiffManifoldBSplineConfigurationDerived, TDiffManifoldBSplineConfigurationDerived> Conf;
			typedef ::bsplines::DiffManifoldBSpline<Conf> BSpline;

			typedef typename TDiffManifoldBSplineConfigurationDerived::ManifoldConf ManifoldConf;

			DesignVariableSegmentBSplineConf(TDiffManifoldBSplineConfigurationDerived derivedConf = TDiffManifoldBSplineConfigurationDerived()) : ParentConf(derivedConf) {}
		};
	}
}

namespace bsplines {
namespace internal {
	template <typename TDiffManifoldBSplineConfiguration, typename TDiffManifoldBSplineConfigurationDerived>
	struct SegmentData< ::aslam::splines::DesignVariableSegmentBSplineConf<TDiffManifoldBSplineConfiguration, TDiffManifoldBSplineConfigurationDerived> > : public SegmentData<TDiffManifoldBSplineConfigurationDerived>, aslam::backend::DesignVariable{
	private:
		typedef SegmentData<TDiffManifoldBSplineConfigurationDerived> parent_t;
		typedef SegmentData< ::aslam::splines::DesignVariableSegmentBSplineConf<TDiffManifoldBSplineConfiguration, TDiffManifoldBSplineConfigurationDerived> > this_t;

	public:
		typedef typename TDiffManifoldBSplineConfiguration::Manifold Manifold;
		typedef typename TDiffManifoldBSplineConfiguration::Manifold::point_t point_t;
		typedef typename TDiffManifoldBSplineConfiguration::Manifold::tangent_vector_t tangent_vector_t;

	private:
		point_t _p_v;
		typename TDiffManifoldBSplineConfiguration::Dimension _dimension;
		typedef typename manifolds::internal::DiffManifoldPointUpdateTraits<typename Manifold::configuration_t> UpdateTraits;
		const Manifold & _manifold;

		enum {
			Dimension = TDiffManifoldBSplineConfiguration::Dimension::VALUE
		};

	public:
		SegmentData(const TDiffManifoldBSplineConfiguration & configuration, const Manifold & manifold, const time_t & time, const point_t &point) : parent_t(configuration, manifold, time, point), _dimension(configuration.getDimension()), _manifold(manifold) {}

		virtual int minimalDimensionsImplementation() const { return _dimension.getValue(); };

		/// \brief Update the design variable.
		virtual void updateImplementation(const double * dp, int size){
			if(!TDiffManifoldBSplineConfiguration::Dimension::IS_DYNAMIC){
				SM_ASSERT_EQ_DBG(std::runtime_error, Dimension, size, "");
			}
			_p_v = this->getControlVertex();
			Eigen::Map<const tangent_vector_t> dpV(dp, size);
			UpdateTraits::update(_manifold, this->getControlVertex(), dpV);
		};

		/// \brief Revert the last state update.
		virtual void revertUpdateImplementation() { this->getControlVertex() = _p_v ; }

		/// Returns the content of the design variable
		virtual void getParametersImplementation(Eigen::MatrixXd& value) const {
			value = this->getControlVertex();
		}

		/// Sets the content of the design variable
		virtual void setParametersImplementation(const Eigen::MatrixXd& value) {
			_p_v = this->getControlVertex();
			this->getControlVertex() = value;
		}

		inline aslam::backend::DesignVariable & getDesignVariable(){ return *this; }
		inline const aslam::backend::DesignVariable & getDesignVariable() const { return *this; }

		/// Computes the minimal distance in tangent space between the current value of the DV and xHat
		virtual void minimalDifferenceImplementation(const Eigen::MatrixXd& xHat, Eigen::VectorXd& outDifference) const;

		/// Computes the minimal distance in tangent space between the current value of the DV and xHat and the jacobian
		virtual void minimalDifferenceAndJacobianImplementation(const Eigen::MatrixXd& xHat, Eigen::VectorXd& outDifference, Eigen::MatrixXd& outJacobian) const;
	};

	template <typename TimePolicy, bool IsInteger = std::numeric_limits<typename TimePolicy::time_t>::is_integer>
	struct TimeExpressionTraits {
		typedef aslam::backend::GenericScalarExpression<typename TimePolicy::time_t> type;
	};

	template <typename TimePolicy>
	struct TimeExpressionTraits<TimePolicy, true> {
		typedef aslam::backend::GenericScalarExpression<aslam::backend::FixedPointNumber<typename TimePolicy::time_t, TimePolicy::getOne()>> type;
	};
}

template <typename TModifiedConf, typename TModifiedDerivedConf>
class DiffManifoldBSpline<aslam::splines::DesignVariableSegmentBSplineConf<TModifiedConf, TModifiedDerivedConf>, aslam::splines::DesignVariableSegmentBSplineConf<TModifiedDerivedConf> > : public bsplines::DiffManifoldBSpline<TModifiedDerivedConf, aslam::splines::DesignVariableSegmentBSplineConf<TModifiedDerivedConf> > {
 public:
	typedef aslam::splines::DesignVariableSegmentBSplineConf<TModifiedDerivedConf> configuration_t;
	typedef bsplines::DiffManifoldBSpline<TModifiedDerivedConf, configuration_t> parent_t;
	typedef typename configuration_t::BSpline spline_t;
	typedef typename parent_t::TimePolicy TimePolicy;
	typedef typename parent_t::time_t time_t;
	typedef typename parent_t::point_t point_t;
	typedef typename parent_t::SegmentIterator SegmentIterator;
	typedef typename parent_t::SegmentConstIterator SegmentConstIterator;
	typedef typename internal::TimeExpressionTraits<TimePolicy>::type TimeExpression;

	typedef spline_t TYPE;

	enum {
		PointSize = TModifiedDerivedConf::ManifoldConf::PointSize::VALUE,
		Dimension = TModifiedDerivedConf::ManifoldConf::Dimension::VALUE
	};

	typedef aslam::backend::DesignVariable dv_t;
	typedef aslam::backend::VectorExpression<PointSize> expression_t;


	DiffManifoldBSpline(const configuration_t & config = configuration_t()) : parent_t(config){}

	size_t numDesignVariables();
	dv_t * designVariable(size_t i);

	void init();

	const std::vector<dv_t *> & getDesignVariables();
	std::vector<dv_t *> getDesignVariables(time_t time);

	// add one Segment at the end of the Spline
	time_t appendSegments(KnotGenerator<time_t> & knotGenerator, int numSegments, const point_t * value);
	void removeSegment();

	template <typename FactoryData_>
	class ExpressionFactory {
	 public:
		typedef boost::shared_ptr<FactoryData_> DataSharedPtr;
		inline const typename FactoryData_::eval_t & getEvaluator() const { return _dataPtr->getEvaluator(); };
		/// \brief get an value expression
		expression_t getValueExpression(int derivativeOrder = 0) const;
	 protected:
		inline const DataSharedPtr & getDataPtr() const { return _dataPtr; }
		inline ExpressionFactory(const FactoryData_ & factoryData) : _dataPtr(new FactoryData_(factoryData)) {}
		friend class DiffManifoldBSpline;
	 private:
		DataSharedPtr _dataPtr;
	};

 protected:
	template<int IMaxDerivativeOrder>
	class ConstTimeFactoryData {
	 public:
		typedef typename spline_t::template Evaluator<IMaxDerivativeOrder> eval_t;
		inline const eval_t & getEvaluator() const { return _eval; };
		inline const spline_t & getSpline() const { return _eval.getSpline(); };
		inline ConstTimeFactoryData(const spline_t & spline, time_t t) : _eval(spline, t) {}
		inline void getDesignVariables(aslam::backend::DesignVariable::set_t & designVariables) const { for(auto & i : _eval) { designVariables.insert(const_cast<aslam::backend::DesignVariable *>(&i.getDesignVariable())); }}
		inline bool hasTimeExpression(){ return false; }
		inline const TimeExpression & getTimeExpression(){ throw std::runtime_error("not implemented"); }
	 protected:
		inline ConstTimeFactoryData() = default;
	 private:
		eval_t _eval;
	};

	template<int IMaxDerivativeOrder>
	class TimeExpressionFactoryData {
	 public:
		typedef typename spline_t::template Evaluator<IMaxDerivativeOrder != Eigen::Dynamic ? IMaxDerivativeOrder + 1 : Eigen::Dynamic> eval_t;
		inline TimeExpressionFactoryData(const spline_t & spline, const TimeExpression & timeExp, time_t lowerBound, time_t upperBound) : _timeExp(timeExp), _lowerBound(lowerBound), _upperBound(upperBound), evalPtr(nullptr), _spline(spline) { SM_ASSERT_LT(std::runtime_error, lowerBound, upperBound, "upper bound must be greater than lower bound."); }
		~TimeExpressionFactoryData(){ delete evalPtr; }
		inline const spline_t & getSpline() const { return _spline; };
		inline void getDesignVariables(aslam::backend::DesignVariable::set_t & designVariables) const;
		const eval_t & getEvaluator() const;
		inline bool hasTimeExpression(){ return true; }
		inline const TimeExpression & getTimeExpression(){ return _timeExp; }
		inline SegmentConstIterator begin() const { return _spline.getFirstRelevantSegmentByLast(_spline.getSegmentIterator(std::max(_spline.getMinTime(), _lowerBound))); }
		inline SegmentConstIterator end() const { auto end = _spline.getSegmentIterator(std::min(_spline.getMaxTime(),_upperBound)); if(end != _spline.end()) end++; return end; }
	 private:
		TimeExpression _timeExp;
		time_t _lowerBound;
		time_t _upperBound;
		mutable eval_t * evalPtr;
		const spline_t & _spline;
	};

 public:
	template <int IMaxDerivativeOrder> inline ExpressionFactory<ConstTimeFactoryData<IMaxDerivativeOrder> > getExpressionFactoryAt(const time_t & t) const {
		return ExpressionFactory<ConstTimeFactoryData<IMaxDerivativeOrder> >(ConstTimeFactoryData<IMaxDerivativeOrder>(*this, t));
	}
	template <int IMaxDerivativeOrder> inline ExpressionFactory<TimeExpressionFactoryData<IMaxDerivativeOrder> > getExpressionFactoryAt(const TimeExpression & t, time_t lowerBound, time_t upperBound) const {
		return ExpressionFactory<TimeExpressionFactoryData<IMaxDerivativeOrder> >(TimeExpressionFactoryData<IMaxDerivativeOrder>(*this, t, lowerBound, upperBound));
	}

 private:
	/// \brief the vector of design variables.
	std::vector< dv_t * > _designVariables;
};
}

namespace aslam {
namespace splines {
	template <typename TBSplineConf>
	class OPTBSpline {
	public:
		typedef aslam::splines::DesignVariableSegmentBSplineConf<TBSplineConf> CONF;
//		typedef ::bsplines::DiffManifoldBSpline<CONF> BSpline;
		typedef typename CONF::BSpline BSpline;
	};
} // namespace splines
} // namespace aslam



#include "implementation/OPTBSplineImpl.hpp"

#endif /* OPTBSPLINE_HPP_ */
