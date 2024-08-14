// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <raptor/upgrade/upgrade.hpp>
#include <raptor/upgrade/upgrade_index.hpp>
#include <raptor/upgrade/load_hibf.hpp>

namespace raptor
{

void raptor_upgrade(upgrade_arguments const & arguments)
{    if (arguments.is_hibf)
    { std::cout<<"test";
        if (arguments.compressed)
            load_hibf<true>(arguments);
        else
            load_hibf<false>(arguments);
    }else{
        if (arguments.compressed)
            upgrade_index<true>(arguments);
        else
            upgrade_index<false>(arguments);
    }
    return;
}

} // namespace raptor
