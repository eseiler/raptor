// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::raptor_search.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <raptor/search/search_hibf.hpp>
#include <raptor/search/search_ibf.hpp>
#include <raptor/search/search_partitioned_ibf.hpp>

void raptor_search_fpga_oneapi(raptor::search_fpga_arguments const & arguments);

namespace raptor
{

void raptor_search(search_arguments const & arguments)
{
    arguments.complete_search_timer.start();

    raptor_search_fpga_oneapi(raptor::search_fpga_arguments{arguments});

    // if (arguments.is_hibf)
    //     search_hibf(arguments);
    // else if (arguments.parts == 1u)
    //     search_ibf(arguments);
    // else
    //     search_partitioned_ibf(arguments);

    arguments.complete_search_timer.stop();

    return;
}

} // namespace raptor
