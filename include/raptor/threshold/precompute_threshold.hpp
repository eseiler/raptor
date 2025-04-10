// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides raptor::threshold::precompute_threshold.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <raptor/threshold/threshold_parameters.hpp>

namespace raptor::threshold
{

[[nodiscard]] std::vector<size_t> precompute_threshold(threshold_parameters const & arguments);

} // namespace raptor::threshold
