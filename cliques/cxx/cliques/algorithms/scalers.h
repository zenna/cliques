/* Copyright (c) Zenna Tavares - zennatavares@gmail.com, 2010-2011 */
#ifndef CLIQUES_SCALERS_H
#define CLIQUES_SCALERS_H

#include <vector>
#include <limits>

namespace cliques {
/**
@brief  A linear scaler functor

This scaler takes (or generates) a range to operate within and scales any new
values to a within a new range
So for example if you have weights of ints from 0 - 344, and you want a float
between 0 and 1, you can use this functor
*/
class linear_scaler {
private:
	float min_value, max_value, shift, range;
	bool do_update_ranges;

public:
	linear_scaler(float min, float max) : min_value(min), max_value(max), do_update_ranges(true){
		update_ranges();
	}

	linear_scaler() {
		min_value = std::numeric_limits<float>::max();
		max_value = -std::numeric_limits<float>::max();
		update_ranges();
	}

	/**
	@brief Use this to update min_value and max_value based on data

	Typical use case is to iterate through your data, calling this function
	with every value as an argument.
	*/
	void update_extrema(float update_value) {
		update_maxima(update_value);
		update_minima(update_value);
	}

	void update_maxima(float update_value) {
		if (update_value > max_value) {
			max_value = update_value;
			do_update_ranges = true;
		}
	}

	void update_minima(float update_value) {
		if (update_value < min_value) {
			min_value = update_value;
			do_update_ranges = true;
		}
	}

	void update_ranges() {
		shift = 0.0 - min_value;
		range = max_value - min_value;
		do_update_ranges = false;
	}

	float operator () (float value) {
		if (do_update_ranges == true) {
			update_ranges();
		}
		return (value + shift)/range;
	}
};

}

#endif // CLIQUES_HSG_H
