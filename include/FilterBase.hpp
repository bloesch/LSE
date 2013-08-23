/*!
* @file 	FilterBase.hpp
* @author 	Michael Blösch
* @date		10.10.2012
* @brief	VUKF Filter for legged robots
 */

#ifndef FILTERBASE_HPP_
#define FILTERBASE_HPP_

#include "Common.hpp"
#include <Eigen/Dense>
#include "Rotations.hpp"

namespace LSE {

/*! Observability Constrained Extended Kalman Filter */
class FilterBase{
public:
	/* -------------------- Constructor/Destructor --------------------- */
	/*! Constructor */
	FilterBase(){};
	/*! Destructor */
	virtual ~FilterBase(){};

	/* -------------------- Filter handling --------------------- */
	/*! Updates the filter to time t
	 * @param[in]	t	desired update time
	 */
	virtual void update(const double& t) = 0;
	/*! Updates the filter to the newest measurement time */
	virtual void update() = 0;
	/*! Return current estimate of robot state (main body)
	 * @return	current robot state
	 */
	virtual State getEst() = 0;
	/*! Return slippage detection
	 * @return	current slippage detection
	 */
	virtual SlippageDetection getSlippage() = 0;
	/*! Resets the filter
	 * @param[in]	t	time used to initialize new state estimate
	 */
	virtual void resetEstimate(const double& t) = 0;
	/*! Returns a string describing the main filter parameters
	 * @param[out] str	string characterize the parameter set of the filter
	 */
	virtual std::string getKeyString() = 0;
};

}

#endif /* FILTERBASE_HPP_ */
