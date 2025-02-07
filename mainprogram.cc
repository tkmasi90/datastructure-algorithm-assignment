// The main program (provided by the course), COMP.CS.300
//
// DO ****NOT**** EDIT THIS FILE!
// (Preferably do not edit this even temporarily. And if you still decide to do so
//  (for debugging, for example), DO NOT commit any changes to git, but revert all
//  changes later. Otherwise you won't be able to get any updates/fixes to this
//  file from git!)

#include <string>
using std::string;
using std::getline;

#include <iostream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::flush;
using std::noskipws;

#include <iomanip>
using std::setfill;
using std::setw;

#include <istream>
using std::istream;

#include <ostream>
using std::ostream;

#include <fstream>
using std::ifstream;

#include <sstream>
using std::istringstream;
using std::ostringstream;
using std::stringstream;

#include <iomanip>
using std::setw;

#include <tuple>
using std::tuple;
using std::make_tuple;
using std::get;
using std::tie;

#include <regex>
using std::regex_match;
using std::regex_search;
using std::smatch;
using std::regex;
using std::sregex_token_iterator;

#include <algorithm>
using std::find_if;
using std::find;
using std::binary_search;
using std::max_element;
using std::max;
using std::min;
using std::shuffle;
using std::sort;

#include <random>
using std::minstd_rand;
using std::uniform_int_distribution;

#include <chrono>

#include <functional>
using std::function;
using std::equal_to;

#include <vector>
using std::vector;

#include <set>
using std::set;

#include <array>
using std::array;

#include <bitset>
using std::bitset;

#include <iterator>
using std::next;

#include <ctime>
using std::time;

#include <memory>
using std::move;

#include <utility>
using std::pair;
using std::make_pair;

#include <cmath>
using std::abs;

#include <cstdlib>
using std::div;

#include <algorithm>
using std::transform;

#include <iterator>
using std::back_inserter;

#include <cstddef>
#include <cassert>


#include "mainprogram.hh"

#include "datastructures.hh"

#ifdef GRAPHICAL_GUI
#include "mainwindow.hh"
#endif

string const MainProgram::PROMPT = "> ";

void MainProgram::test_get_functions(AffiliationID id)
{
    ds_.get_affiliation_name(id);
    ds_.get_affiliation_coord(id);
}

MainProgram::CmdResult MainProgram::cmd_add_affiliation(ostream& /*output*/, MatchIter begin, MatchIter end)
{
    AffiliationID id = *begin++;
    string name = *begin++;
    string xstr = *begin++;
    string ystr = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    int x = convert_string_to<int>(xstr);
    int y = convert_string_to<int>(ystr);

    bool success = ds_.add_affiliation(id, name, {x, y});

    view_dirty = true;
    return {ResultType::IDLIST, CmdResultIDs{{}, {success ? id : NO_AFFILIATION}}};
}

MainProgram::CmdResult MainProgram::cmd_affiliation_info(std::ostream& /*output*/, MatchIter begin, MatchIter end)
{
    AffiliationID id = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    return {ResultType::IDLIST, CmdResultIDs{{}, {id}}};
}

void MainProgram::test_affiliation_info()
{
    if (random_affiliations_added_ > 0) // Don't do anything if there's no affiliations
    {
        auto id = random_affiliation();
        test_get_functions(id);
    }
}

MainProgram::CmdResult MainProgram::cmd_change_affiliation_coord(std::ostream& /*output*/, MainProgram::MatchIter begin, MainProgram::MatchIter end)
{
    AffiliationID id = *begin++;
    string xstr = *begin++;
    string ystr = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    int x = convert_string_to<int>(xstr);
    int y = convert_string_to<int>(ystr);

    bool success = ds_.change_affiliation_coord(id, {x,y});

    view_dirty = true;
    return {ResultType::IDLIST, CmdResultIDs{{}, {success ? id : NO_AFFILIATION}}};
}

MainProgram::CmdResult MainProgram::cmd_get_publications_after(std::ostream &output, MatchIter begin, MatchIter end)
{
    AffiliationID affiliationid = *begin++;
    Year time = convert_string_to<Year>(*begin++);
    assert( begin == end && "Impossible number of parameters!");

    auto publications = ds_.get_publications_after(affiliationid, time);
    
    if (publications.size() == 1 && publications.front() == std::make_pair(NO_YEAR, NO_PUBLICATION))
    {
        output << "No such publications found (NO_YEAR, NO_PUBLICATION returned)" << endl;
        return {};
    }

    if (!publications.empty())
    {
        output << "Publications from affiliation ";
        print_affiliation_brief(affiliationid, output, false);
        output << " after year " << setw(4) << setfill('0') << time << ":" << endl;
        for (auto& [deptime, publicationid] : publications)
        {
            output << " " << publicationid << " at " << setw(4) << setfill('0') << deptime << endl;
        }
    }
    else
    {
        output << "No publications from affiliation ";
        print_affiliation_brief(affiliationid, output, false);
        output << " after year " << time << endl;
    }

    return {};
}

void MainProgram::test_get_publications_after()
{
    if (random_affiliations_added_ > 0) // Don't do anything if there's no affiliations
    {
        auto id = random_affiliation();
        ds_.get_publications_after(id, get_random_year());
    }
}

void MainProgram::test_change_affiliation_coord()
{
    if (random_affiliations_added_ > 0) // Don't do anything if there's no affiliations
    {
        auto id = random_affiliation();
        ds_.change_affiliation_coord(id, get_random_coords());
    }
}

MainProgram::CmdResult MainProgram::cmd_add_reference(std::ostream& output, MainProgram::MatchIter begin, MainProgram::MatchIter end)
{
    // TODO check order of parameters!!
    PublicationID id = convert_string_to<PublicationID>(*begin++);
    PublicationID parentid = convert_string_to<PublicationID>(*begin++);
    assert( begin == end && "Impossible number of parameters!");

    bool ok = ds_.add_reference(id, parentid);
    if (ok)
    {
        try
        {
            auto referencename = ds_.get_publication_name(id);
            auto parentname = ds_.get_publication_name(parentid);
            output << "Added '" << referencename << "' as a reference of '" << parentname << "'" << endl;
        }
        catch (NotImplemented&)
        {
            output << "Added a reference to a publication." << endl;
        }
        return {ResultType::IDLIST, CmdResultIDs{{id, parentid}, {}}};
    }
    else
    {
        output << "Adding a reference failed!" << endl;
        return {};
    }
}

MainProgram::CmdResult MainProgram::cmd_add_affiliation_to_publication(std::ostream &output, MatchIter begin, MatchIter end)
{
    AffiliationID affiliationid = *begin++;
    PublicationID publicationid = convert_string_to<PublicationID>(*begin++);
    assert( begin == end && "Impossible number of parameters!");

    bool ok = ds_.add_affiliation_to_publication(affiliationid, publicationid);
    if (ok)
    {
        try
        {
            auto affiliationname = ds_.get_affiliation_name(affiliationid);
            auto publicationname = ds_.get_publication_name(publicationid);
            output << "Added '" << affiliationname << "' as an affiliation to publication '" << publicationname << "'" << endl;
        }
        catch (NotImplemented&)
        {
            output << "Added a affiliation to publication." << endl;
        }
        return {ResultType::IDLIST, CmdResultIDs{{publicationid}, {affiliationid}}};
    }
    else
    {
        output << "Adding a affiliation to publication failed!" << endl;
        return {};
    }
}

MainProgram::CmdResult MainProgram::cmd_get_affiliations_closest_to(std::ostream &output, MatchIter begin, MatchIter end)
{
    string xstr = *begin++;
    string ystr = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    int x = convert_string_to<int>(xstr);
    int y = convert_string_to<int>(ystr);

    auto affiliations = ds_.get_affiliations_closest_to({x,y});
    if (affiliations.empty())
    {
        output << "No affiliations!" << endl;
    }

    return {ResultType::IDLIST, CmdResultIDs{{}, affiliations}};
}

MainProgram::CmdResult MainProgram::cmd_get_closest_common_parent(std::ostream &output, MatchIter begin, MatchIter end)
{
    PublicationID publicationid1 = convert_string_to<PublicationID>(*begin++);
    PublicationID publicationid2 = convert_string_to<PublicationID>(*begin++);
    assert( begin == end && "Impossible number of parameters!");


    auto publicationid = ds_.get_closest_common_parent(publicationid1, publicationid2);
    if (publicationid == NO_PUBLICATION)
    {
        output << "No common referring publication found." << endl;
    }

    return {ResultType::IDLIST, CmdResultIDs{{publicationid1, publicationid2, publicationid}, {}}};
}

MainProgram::CmdResult MainProgram::cmd_remove_publication(std::ostream &output, MatchIter begin, MatchIter end)
{
    PublicationID pubid = convert_string_to<PublicationID>(*begin++);
    assert( begin == end && "Impossible number of parameters!");

    auto pubname = ds_.get_publication_name(pubid);
    bool success = ds_.remove_publication(pubid);
    if (success)
    {
        output << pubname << " removed." << endl;
        view_dirty = true;
        return {};
    }
    else
    {
        return {ResultType::IDLIST, CmdResultIDs{{NO_PUBLICATION}, {}}};
    }
}

MainProgram::CmdResult MainProgram::cmd_get_parent(std::ostream& output, MatchIter begin, MatchIter end)
{
    PublicationID pubid = convert_string_to<PublicationID>(*begin++);
    assert( begin == end && "Impossible number of parameters!");

    auto refid = ds_.get_parent(pubid);
    if (refid != NO_PUBLICATION) {
        return {ResultType::IDLIST, CmdResultIDs{{refid},{}}};
    }
    output << "No references or publication doesn't exist." << std::endl;
    return {};
}

MainProgram::CmdResult MainProgram::cmd_get_referenced_by_chain(std::ostream &output, MatchIter begin, MatchIter end)
{
    PublicationID pubid = convert_string_to<PublicationID>(*begin++);
    assert( begin == end && "Impossible number of parameters!");
    auto reference_chain = ds_.get_referenced_by_chain(pubid);
    if (reference_chain.empty()) { output << "Publication is not cited anywhere." << std::endl; }
    return {ResultType::IDLIST, CmdResultIDs{reference_chain, {}}};
}

MainProgram::CmdResult MainProgram::cmd_get_direct_references(std::ostream &output, MatchIter begin, MatchIter end)
{
    PublicationID pubid = convert_string_to<PublicationID>(*begin++);
    assert( begin == end && "Impossible number of parameters!");
    auto references = ds_.get_direct_references(pubid);
    if (references.empty()) { output << "Publication has no direct references." << std::endl; }
    std::sort(references.begin(),references.end());
    return {ResultType::IDLIST, CmdResultIDs{references, {}}};
}

MainProgram::CmdResult MainProgram::cmd_get_publications(std::ostream& output, MainProgram::MatchIter begin, MainProgram::MatchIter end)
{
    AffiliationID id = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    auto result = ds_.get_publications(id);
    if (result.empty()) { output << "Affiliation has no publications." << std::endl; }
    std::sort(result.begin(),result.end());
    return {ResultType::IDLIST, CmdResultIDs{result, {id}}};
}

void MainProgram::test_get_publications()
{
    if (random_affiliations_added_ > 0) // Don't do anything if there's no affiliations
    {
        auto id = random_affiliation();
        ds_.get_publications(id);
    }
}

void MainProgram::test_get_all_references()
{
    if (random_publications_added_ > 0) // Don't do anything if there's no publications
    {
        auto id = random_root_publication();
        ds_.get_all_references(id);
    }
}

MainProgram::CmdResult MainProgram::cmd_remove_affiliation(ostream& output, MatchIter begin, MatchIter end)
{
    string id = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    auto name = ds_.get_affiliation_name(id);
    bool success = ds_.remove_affiliation(id);
    if (success)
    {
        output << name << " removed." << endl;
        view_dirty = true;
        return {};
    }
    else
    {
        return {ResultType::IDLIST, CmdResultIDs{{}, {NO_AFFILIATION}}};
    }
}

void MainProgram::test_remove_affiliation()
{
    // Choose random number to remove
    if (random_affiliations_added_ > 0) // Don't remove if there's nothing to remove
    {
        auto affiliationid = random_affiliation(); // there is a risk of getting the same id to remove multiple times -> usually takes less time than with existing affiliation
        ds_.remove_affiliation(affiliationid);
    }
}

void MainProgram::add_random_affiliations_publications(unsigned int size, Coord min, Coord max, const std::vector<Coord> &coordinates)
{
    if(coordinates.size()!=size){
        for (unsigned int i = 0; i < size; ++i)
        {
            auto name = n_to_name(random_affiliations_added_);
            AffiliationID id = n_to_affiliationid(random_affiliations_added_);

            ds_.add_affiliation(id, name, get_random_coords(min, max));

            ++random_affiliations_added_;
        }
    } else {
        for (unsigned int i = 0; i < size; ++i)
        {
            auto name = n_to_name(random_affiliations_added_);
            AffiliationID id = n_to_affiliationid(random_affiliations_added_);

            ds_.add_affiliation(id, name, coordinates.at(i));

            ++random_affiliations_added_;
        }
    }


    for (unsigned int i = 0; i< size; ++i) {
        auto publicationid = n_to_publicationid(random_publications_added_);

        vector<AffiliationID> affiliations;
        for (int j=0; j<4; ++j)
        {
            affiliations.push_back(random_affiliation());
        }
        ds_.add_publication(publicationid, convert_to_string(publicationid), get_random_year(), std::move(affiliations));

        // Add area as subarea so that we get a binary tree
        if (random_publications_added_ > 0)
        {
            auto parentid = n_to_publicationid(random_publications_added_ / 2);
            ds_.add_reference(publicationid, parentid);
        }
        ++random_publications_added_;
    }
}

MainProgram::CmdResult MainProgram::cmd_random_affiliations(ostream& output, MatchIter begin, MatchIter end)
{
    string sizestr = *begin++;
    string minxstr = *begin++;
    string minystr = *begin++;
    string maxxstr = *begin++;
    string maxystr = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    unsigned int size = convert_string_to<unsigned int>(sizestr);

    Coord min = RANDOM_MIN_COORD;
    Coord max = RANDOM_MAX_COORD;
    std::unordered_set<Coord,CoordHash> exclude_list;
    auto affiliations = ds_.get_all_affiliations();
    for(const auto& affid : affiliations){
        exclude_list.insert(ds_.get_affiliation_coord(affid));
    }
    if (!minxstr.empty() && !minystr.empty() && !maxxstr.empty() && !maxystr.empty())
    {
        min.x = convert_string_to<unsigned int>(minxstr);
        min.y = convert_string_to<unsigned int>(minystr);
        max.x = convert_string_to<unsigned int>(maxxstr);
        max.y = convert_string_to<unsigned int>(maxystr);
    }
    else
    {
        if (!exclude_list.empty())
        {
            // Find out the min and max coordinates current affiliations reside within -> generates more data within the same region
            min = {std::numeric_limits<int>::max(), std::numeric_limits<int>::max()};
            max = {std::numeric_limits<int>::min(), std::numeric_limits<int>::min()};
            for (auto const& coord : exclude_list)
            {
                auto x = coord.x;
                auto y = coord.y;
                if (x < min.x) { min.x = x; }
                if (y < min.y) { min.y = y; }
                if (x > max.x) { max.x = x; }
                if (y > max.y) { max.y = y; }
            }
        }
    }

    std::vector<Coord> unique_new_coords;
    try {
        unique_new_coords = get_unique_coords(size,exclude_list,min,max);
    } catch (...) {
        output << "Impossible to create such number of unique coordinates within perimeters" <<endl;
        return{};
    }
    add_random_affiliations_publications(size, min, max, unique_new_coords);

    output << "Added: " << size << " affiliations and publications." << endl;

    view_dirty = true;

    return {};
}

void MainProgram::test_random_affiliations()
{
    add_random_affiliations_publications(1);
}

void MainProgram::test_remove_publication()
{
    if (random_publications_added_ > 0){
        auto publicationid = random_publication(); // same as with remove_affiliation, could result in already removed id
        ds_.remove_publication(publicationid);
    }
}

void MainProgram::test_get_parent()
{
    if (random_publications_added_ > 0){
        auto publicationid = random_publication();
        ds_.get_parent(publicationid);
    }
}

void MainProgram::test_get_referenced_by_chain()
{
    if (random_publications_added_ > 0){
        auto publicationid = random_leaf_publication();
        ds_.get_referenced_by_chain(publicationid);
    }
}

void MainProgram::test_get_direct_references()
{
    if (random_publications_added_ > 0) {
        auto publicationid = random_publication();
        ds_.get_direct_references(publicationid);
    }
}

void MainProgram::test_get_affiliations()
{
    if (random_publications_added_ > 0) {
        auto publicationid = random_publication();
        ds_.get_affiliations(publicationid);
    }
}

void MainProgram::test_get_affiliation_count()
{
    ds_.get_affiliation_count();
}

void MainProgram::test_get_all_publications()
{
    ds_.all_publications();
}

void MainProgram::test_add_affiliation_to_publication()
{
    if (random_publications_added_ > 0 || random_affiliations_added_ > 0) {
        auto publicationid = random_publication();
        auto affiliationid = random_affiliation();
        ds_.add_affiliation_to_publication(affiliationid, publicationid);
    }
}

void MainProgram::test_get_connected_affiliations()
{
    if (random_publications_added_ > 0 ){
        auto affiliationid = random_affiliation();
        ds_.get_connected_affiliations(affiliationid);
    }
}

void MainProgram::test_get_all_connections()
{
    ds_.get_all_connections();
}

void MainProgram::test_get_any_path()
{
    if (random_publications_added_ > 0 ){
        auto fromid = random_affiliation();
        auto toid = random_affiliation();
        ds_.get_any_path(fromid, toid);
    }
}

void MainProgram::test_get_path_with_least_affiliations()
{
    if (random_publications_added_ > 0 ){
        auto fromid = random_affiliation();
        auto toid = random_affiliation();
        ds_.get_path_with_least_affiliations(fromid, toid);
    }
}

void MainProgram::test_get_path_of_least_friction()
{
    if (random_publications_added_ > 0 ){
        auto fromid = random_affiliation();
        auto toid = random_affiliation();
        ds_.get_path_of_least_friction(fromid, toid);
    }
}

void MainProgram::test_get_shortest_path()
{
    if (random_publications_added_ > 0 ){
        auto fromid = random_affiliation();
        auto toid = random_affiliation();
        ds_.get_shortest_path(fromid, toid);
    }
}

Coord MainProgram::get_random_coords(const Coord min, const Coord max)
{
    int x = random<int>(min.x, max.x);
    int y = random<int>(min.y, max.y);
    return {x, y};
}

Year MainProgram::get_random_year(const Year min, const Year max)
{
    return random<int>(min, max);
}

std::vector<Coord> MainProgram::get_unique_coords(const unsigned int n, const std::unordered_set<Coord, CoordHash> &exclude_list, const Coord min, const Coord max)
{
    const unsigned int max_unique_coords=abs(max.x-min.x)*abs(max.y-min.y)-exclude_list.size();
    if(n>max_unique_coords){
        throw NotImplemented("Impossible to create such number of unique coordinates within perimeters");
    }
    //assert(n<max_unique_coords&&"Impossible to create such number of unique coordinates within perimeters");

    if(n>max_unique_coords/2){
        // let's not use the lottery in this case, let's list the coordinates not in use
        std::vector<Coord> retvec;
        retvec.reserve(max_unique_coords);
        for(auto x = min.x; x<max.x;++x){
            for(auto y=min.y;y<max.y;++y){
                Coord newCoord ={x,y};
                if(exclude_list.find(newCoord)==exclude_list.end()){
                    retvec.push_back(newCoord);
                }
            }
        }
        // and get the required number of random_coordinates
        std::shuffle(retvec.begin(),retvec.end(),rand_engine_);
        retvec.erase(std::next(retvec.begin(),n),retvec.end());
        return retvec;
    } else {
        std::unordered_set<Coord,CoordHash> coords;

        coords.reserve(n);
        while(coords.size()<n){
            Coord newCoord = get_random_coords(min, max);
            if(exclude_list.find(newCoord)==exclude_list.end()){
                coords.insert(newCoord);
            }
        }
        std::vector<Coord> retvec(coords.begin(),coords.end());
        return retvec;
    }


}

MainProgram::CmdResult MainProgram::cmd_get_affiliation_count(ostream& output, MatchIter begin, MatchIter end)
{
    assert( begin == end && "Impossible number of parameters!");

    output << "Number of affiliations: " << ds_.get_affiliation_count() << endl;

    return {};
}

MainProgram::CmdResult MainProgram::cmd_get_all_affiliations(ostream& output, MatchIter begin, MatchIter end)
{
    assert( begin == end && "Impossible number of parameters!");

    auto affiliations = ds_.get_all_affiliations();
    if (affiliations.empty())
    {
        output << "No affiliations!" << endl;
    }

    std::sort(affiliations.begin(), affiliations.end());
    return {ResultType::IDLIST, CmdResultIDs{{}, affiliations}};
}

MainProgram::CmdResult MainProgram::cmd_add_publication(std::ostream& /*output*/, MatchIter begin, MatchIter end)
{
    PublicationID id = convert_string_to<PublicationID>(*begin++);
    string name = *begin++;
    Year year = convert_string_to<Year>(*begin++);
    string affilsstr = *begin++;

    assert( begin == end && "Impossible number of parameters!");

    vector<AffiliationID> affiliations;

    smatch affil;
    auto sbeg = affilsstr.cbegin();
    auto send = affilsstr.cend();
    for ( ; regex_search(sbeg, send, affil, affil_regex_); sbeg = affil.suffix().first)
    {
        affiliations.push_back(affil[1]);
    }
    bool success = ds_.add_publication(id, name, year, affiliations);

    view_dirty = true;
    return {ResultType::IDLIST, CmdResultIDs{{success ? id : NO_PUBLICATION}, {}}};
}

MainProgram::CmdResult MainProgram::cmd_get_all_publications(std::ostream &output, MatchIter begin, MatchIter end)
{
    assert( begin == end && "Impossible number of parameters!");

    auto publications = ds_.all_publications();
    if (publications.empty())
    {
        output << "No publications!" << endl;
    }

    std::sort(publications.begin(), publications.end());
    return {ResultType::IDLIST, CmdResultIDs{publications, {}}};
}

MainProgram::CmdResult MainProgram::cmd_publication_info(std::ostream& /*output*/, MatchIter begin, MatchIter end)
{
    PublicationID id = convert_string_to<PublicationID>(*begin++);
    assert( begin == end && "Impossible number of parameters!");

    return {ResultType::IDLIST, CmdResultIDs{{id}, {}}};
}

MainProgram::CmdResult MainProgram::cmd_get_all_references(std::ostream &output, MatchIter begin, MatchIter end)
{
    PublicationID publicationid = convert_string_to<PublicationID>(*begin++);
    assert( begin == end && "Impossible number of parameters!");

    auto references = ds_.get_all_references(publicationid);
    if (references.empty())
    {
        output << "No (direct) references!" << endl;
    }

    std::sort(references.begin(), references.end());
    references.insert(references.begin(), publicationid); // Add parameter as the first publication
    return {ResultType::IDLIST, CmdResultIDs{references, {}}};
}

Distance MainProgram::calc_distance(Coord c1, Coord c2)
{
    if (c1 == NO_COORD || c2 == NO_COORD) { return NO_DISTANCE; }

    long long int deltax = c1.x - c2.x;
    long long int deltay = c1.y - c2.y;
    return static_cast<Distance>(std::sqrt(deltax*deltax + deltay*deltay));
}

MainProgram::CmdResult MainProgram::cmd_clear_all(ostream& output, MatchIter begin, MatchIter end)
{
    assert(begin == end && "Invalid number of parameters");

    ds_.clear_all();
    init_primes();

    output << "Cleared all affiliations and publications" << endl;

    view_dirty = true;

    return {};
}

string MainProgram::print_affiliation(AffiliationID id, ostream& output, bool nl)
{
    try
    {
        if (id != NO_AFFILIATION)
        {
            auto name = ds_.get_affiliation_name(id);
            auto xy = ds_.get_affiliation_coord(id);
            if (!name.empty())
            {
                output << name << ": ";
            }
            else
            {
                output << "*: ";
            }

            output << "pos=";
            print_coord(xy, output, false);
            output << ", id=" << id;
            if (nl) { output << endl; }

            ostringstream retstream;
            retstream << id;
            return retstream.str();
        }
        else
        {
            output << "--NO_AFFILIATION--";
            if (nl) { output << endl; }
            return "";
        }
    }
    catch (NotImplemented const& e)
    {
        output << endl << "NotImplemented while printing affiliation : " << e.what() << endl;
        std::cerr << endl << "NotImplemented while printing affiliation : " << e.what() << endl;
        return "";
    }
}

string MainProgram::print_affiliation_brief(AffiliationID id, std::ostream &output, bool nl)
{
    try
    {
        if (id != NO_AFFILIATION)
        {
            auto name = ds_.get_affiliation_name(id);
            if (!name.empty())
            {
                output << name << " ";
            }
            else
            {
                output << "* ";
            }

            output << "(" << id << ")";
            if (nl) { output << endl; }

            ostringstream retstream;
            retstream << id;
            return retstream.str();
        }
        else
        {
            output << "--NO_AFFILIATION--";
            if (nl) { output << endl; }
            return "";
        }
    }
    catch (NotImplemented const& e)
    {
        output << endl << "NotImplemented while printing affiliation : " << e.what() << endl;
        std::cerr << endl << "NotImplemented while printing affiliation : " << e.what() << endl;
        return "";
    }
}

string MainProgram::print_publication(PublicationID id, std::ostream &output, bool nl)
{
    try
    {
        if (id != NO_PUBLICATION)
        {
            auto name = ds_.get_publication_name(id);
            auto year = ds_.get_publication_year(id);
            if (!name.empty())
            {
                output << name << ": ";
            }
            else
            {
                output << "*: ";
            }
            output << "year=";
            if (year == NO_YEAR) {
                output << "--NO_YEAR--";
            } else {
                output << std::to_string(year);
            }
            output << ", id=" << id;
            if (nl) { output << endl; }

            ostringstream retstream;
            retstream << id;
            return retstream.str();
        }
        else
        {
            output << "--NO_PUBLICATION--";
            if (nl) { output << endl; }
            return "";
        }
    }
    catch (NotImplemented const& e)
    {
        output << endl << "NotImplemented while printing publication : " << e.what() << endl;
        std::cerr << endl << "NotImplemented while printing publication : " << e.what() << endl;
        return "";
    }
}

MainProgram::CmdResult MainProgram::cmd_find_affiliation_with_coord(ostream& /* output */, MatchIter begin, MatchIter end)
{
    string xstr = *begin++;
    string ystr = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    int x = convert_string_to<int>(xstr);
    int y = convert_string_to<int>(ystr);

    auto result = ds_.find_affiliation_with_coord({x,y});

    return {ResultType::IDLIST, CmdResultIDs{{}, {result}}};
}

MainProgram::CmdResult MainProgram::cmd_get_affiliations(std::ostream &output, MatchIter begin, MatchIter end)
{
    auto pubid = convert_string_to<PublicationID>(*begin++);
    assert( begin == end && "Impossible number of parameters!");

    auto affiliations = ds_.get_affiliations(pubid);
    std::sort(affiliations.begin(),affiliations.end());
    if (affiliations.empty()) { output << "Publication has no affiliations." << std::endl; }
    return {ResultType::IDLIST, CmdResultIDs{{pubid},affiliations}};
}

MainProgram::CmdResult MainProgram::cmd_get_connected_affiliations(std::ostream &output, MatchIter begin, MatchIter end)
{
    auto affid = convert_string_to<AffiliationID>(*begin++);
    assert( begin == end && "Impossible number of parameters!");
    auto connections = ds_.get_connected_affiliations(affid);

    if (connections.empty()){
        output << "No connections from " <<affid<<"!"<< endl;
        return {};
    }
    std::sort(connections.begin(),connections.end(),[](auto& c1, auto& c2){return c1.weight == c2.weight ? (c1.aff2 <c2.aff2) : c1.weight > c2.weight;});
    return {ResultType::NEIGHBOURLIST,connections};
}

MainProgram::CmdResult MainProgram::cmd_get_all_connections(std::ostream &output, MatchIter begin, MatchIter end)
{
    assert( begin == end && "Impossible number of parameters!");
    auto connections = ds_.get_all_connections();
    if (connections.empty()){
        output << "No connections!" << endl;
        return {};
    }
    std::sort(connections.begin(),connections.end(),[](auto& c1, auto& c2){return c1.weight != c2.weight ? c1.weight > c2.weight : (c1.aff1 == c2.aff1 ? c1.aff2 < c2.aff2 : c1.aff1<c2.aff1);});
    return {ResultType::CONNECTIONLIST,connections};
}

MainProgram::CmdResult MainProgram::cmd_get_any_path(std::ostream &output, MatchIter begin, MatchIter end)
{
    auto sourceid = convert_string_to<AffiliationID>(*begin++);
    auto targetid = convert_string_to<AffiliationID>(*begin++);
    assert( begin == end && "Impossible number of parameters!");
    auto route = ds_.get_any_path(sourceid, targetid);
    CmdResultRoute path;
    if (route.empty())
    {
        output << "No route found! (empty route returned)" << endl;
    }
    else {
        Distance distance = 0;
        std::transform(route.begin(),route.end(),std::back_inserter(path),[&distance](auto& connection){
            return std::tuple<AffiliationID, Weight, AffiliationID, Distance>(connection.aff1,connection.weight,connection.aff2,++distance);
        });
    }
    return {ResultType::ROUTE, path};
}

MainProgram::CmdResult MainProgram::cmd_get_path_with_least_affiliations(std::ostream &output, MatchIter begin, MatchIter end)
{
    auto sourceid = convert_string_to<AffiliationID>(*begin++);
    auto targetid = convert_string_to<AffiliationID>(*begin++);
    assert( begin == end && "Impossible number of parameters!");
    auto route = ds_.get_path_with_least_affiliations(sourceid, targetid);
    CmdResultRoute path;
    if (route.empty())
    {
        output << "No route found! (empty route returned)" << endl;
    }
    else {
        Distance distance = 0;
        std::transform(route.begin(),route.end(),std::back_inserter(path),[&distance](auto& connection){
            return std::tuple<AffiliationID, Weight, AffiliationID, Distance>(connection.aff1,connection.weight,connection.aff2,++distance);
        });
    }
    return {ResultType::ROUTE, path};
}

MainProgram::CmdResult MainProgram::cmd_get_path_of_least_friction(std::ostream &output, MatchIter begin, MatchIter end)
{
    auto sourceid = convert_string_to<AffiliationID>(*begin++);
    auto targetid = convert_string_to<AffiliationID>(*begin++);
    assert( begin == end && "Impossible number of parameters!");
    auto route = ds_.get_path_of_least_friction(sourceid, targetid);
    CmdResultRoute path;
    if (route.empty())
    {
        output << "No route found! (empty route returned)" << endl;
    }
    else {
        Distance distance = 0;
        std::transform(route.begin(),route.end(),std::back_inserter(path),[&distance](auto& connection){
            return std::tuple<AffiliationID, Weight, AffiliationID, Distance>(connection.aff1,connection.weight,connection.aff2,++distance);
        });
    }
    return {ResultType::ROUTE, path};
}

MainProgram::CmdResult MainProgram::cmd_get_shortest_path(std::ostream &output, MatchIter begin, MatchIter end)
{
    auto sourceid = convert_string_to<AffiliationID>(*begin++);
    auto targetid = convert_string_to<AffiliationID>(*begin++);
    assert( begin == end && "Impossible number of parameters!");
    auto route = ds_.get_shortest_path(sourceid, targetid);
    CmdResultRoute path;
    if (route.empty())
    {
        output << "No route found! (empty route returned)" << endl;
    }
    else {
        std::transform(route.begin(),route.end(),std::back_inserter(path),[](auto& connection){
            return std::tuple<AffiliationID, Weight, AffiliationID, Distance>(connection.first.aff1,connection.first.weight,connection.first.aff2,connection.second);
        });
    }
    return {ResultType::ROUTE, path};
}

AffiliationID MainProgram::random_affiliation()
{
    return n_to_affiliationid(random<decltype(random_affiliations_added_)>(0, random_affiliations_added_));
}

PublicationID MainProgram::random_publication()
{
    return n_to_publicationid(random<decltype(random_publications_added_)>(0, random_publications_added_));
}

PublicationID MainProgram::random_root_publication()
{
    unsigned long end = ROOT_BIAS_MULTIPLIER * random_publications_added_;
    if (end == 0 ) {
        return 0;
    }
    return n_to_publicationid(random<decltype(random_publications_added_)>(0, end));
}

PublicationID MainProgram::random_leaf_publication()
{
    unsigned long start = LEAF_BIAS_MULTIPLIER * random_publications_added_;
    if (start == random_publications_added_) {
        start = 0;
    }
    return n_to_publicationid(random<decltype(random_publications_added_)>(start, random_publications_added_));
}

void MainProgram::test_find_affiliation_with_coord()
{
    ds_.find_affiliation_with_coord(get_random_coords());
}

void MainProgram::test_publication_info()
{
    if (random_publications_added_ > 0) // Don't do anything if there's no publications
    {
        auto id = random_publication();
        ds_.get_publication_name(id);
        ds_.get_publication_year(id);
    }
}

void MainProgram::test_affiliations_closest_to()
{
    ds_.get_affiliations_closest_to(get_random_coords());
}

void MainProgram::test_get_closest_common_parent()
{
    if (random_publications_added_ > 0) // Don't do anything if there's no publications
    {
        auto id1 = random_leaf_publication();
        auto id2 = random_leaf_publication();
        ds_.get_closest_common_parent(id1, id2);
    }
}

MainProgram::CmdResult MainProgram::cmd_randseed(std::ostream& output, MatchIter begin, MatchIter end)
{
    string seedstr = *begin++;
    assert(begin == end && "Invalid number of parameters");

    unsigned long int seed = convert_string_to<unsigned long int>(seedstr);

    rand_engine_.seed(seed);
    init_primes();

    output << "Random seed set to " << seed << endl;

    return {};
}

MainProgram::CmdResult MainProgram::cmd_read(std::ostream& output, MatchIter begin, MatchIter end)
{
    string filename = *begin++;
    string silentstr =  *begin++;
    assert( begin == end && "Impossible number of parameters!");

    bool silent = !silentstr.empty();
    ostream* new_output = &output;

    ostringstream dummystr; // Given as output if "silent" is specified, the output is discarded
    if (silent)
    {
        new_output = &dummystr;
    }

    ifstream input(filename);
    if (input)
    {
        output << "** Commands from '" << filename << "'" << endl;
        command_parser(input, *new_output, PromptStyle::NORMAL);
        if (silent) { output << "...(output discarded in silent mode)..." << endl; }
        output << "** End of commands from '" << filename << "'" << endl;
    }
    else
    {
        output << "Cannot open file '" << filename << "'!" << endl;
    }

    return {};
}


MainProgram::CmdResult MainProgram::cmd_testread(std::ostream& output, MatchIter begin, MatchIter end)
{
    string infilename = *begin++;
    string outfilename = *begin++;
    assert( begin == end && "Impossible number of parameters!");

    ifstream input(infilename);
    if (input)
    {
        ifstream expected_output(outfilename);
        if (output)
        {
            stringstream actual_output;
            command_parser(input, actual_output, PromptStyle::NO_NESTING);

            vector<string> actual_lines;
            while (actual_output)
            {
                string line;
                getline(actual_output, line);
                if (!actual_output) { break; }
                actual_lines.push_back(line);
            }

            vector<string> expected_lines;
            while (expected_output)
            {
                string line;
                getline(expected_output, line);
                if (!expected_output) { break; }
                expected_lines.push_back(line);
            }

            string heading_actual = "Actual output";
            unsigned int actual_max_length = heading_actual.length();
            auto actual_max_iter = max_element(actual_lines.cbegin(), actual_lines.cend(),
                                               [](string s1, string s2){ return s1.length() < s2.length(); });
            if (actual_max_iter != actual_lines.cend())
            {
                actual_max_length = actual_max_iter->length();
            }

            string heading_expected = "Expected output";
            unsigned int expected_max_length = heading_expected.length();
            auto expected_max_iter = std::max_element(expected_lines.cbegin(), expected_lines.cend(),
                                                      [](string s1, string s2){ return s1.length() < s2.length(); });
            if (expected_max_iter != expected_lines.cend())
            {
                expected_max_length = expected_max_iter->length();
            }

            auto pos_actual = actual_lines.cbegin();
            auto pos_expected = expected_lines.cbegin();
            output << "  " << heading_actual << string(actual_max_length - heading_actual.length(), ' ') << " | " << heading_expected << endl;
            output << "--" << string(actual_max_length, '-') << "-|-" << string(expected_max_length, '-') << endl;

            bool lines_ok = true;
            while (pos_expected != expected_lines.cend() || pos_actual != actual_lines.cend())
            {
                if (pos_expected != expected_lines.cend())
                {
                    if (pos_actual != actual_lines.cend())
                    {
                        bool ok = (*pos_expected == *pos_actual);
                        output << (ok ? ' ' : '?') << ' ' << *pos_actual << string(actual_max_length - pos_actual->length(), ' ')
                               << " | " << *pos_expected << endl;
                        lines_ok = lines_ok && ok;
                        ++pos_actual;
                    }
                    else
                    { // Actual output was too short
                        output << "? " << string(actual_max_length, ' ')
                               << " | " << *pos_expected << endl;
                        lines_ok = false;
                    }
                    ++pos_expected;
                }
                else
                { // Actual output was too long
                    output << "? " << *pos_actual << string(actual_max_length - pos_actual->length(), ' ')
                           << " | " << endl;
                    lines_ok = false;
                    ++pos_actual;
                }
            }
            if (lines_ok)
            {
                output << "**No differences in output.**" << endl;
                if (test_status_ == TestStatus::NOT_RUN)
                {
                    test_status_ = TestStatus::NO_DIFFS;
                }
            }
            else
            {
                output << "**Differences found! (Lines beginning with '?')**" << endl;
                test_status_ = TestStatus::DIFFS_FOUND;
            }

        }
        else
        {
            output << "Cannot open file '" << outfilename << "'!" << endl;
        }
    }
    else
    {
        output << "Cannot open file '" << infilename << "'!" << endl;
    }

    return {};
}

MainProgram::CmdResult MainProgram::cmd_stopwatch(std::ostream& output, MatchIter begin, MatchIter end)
{
    string on = *begin++;
    string off = *begin++;
    string next = *begin++;
    assert(begin == end && "Invalid number of parameters");

    if (!on.empty())
    {
        stopwatch_mode = StopwatchMode::ON;
        output << "Stopwatch: on" << endl;
    }
    else if (!off.empty())
    {
        stopwatch_mode = StopwatchMode::OFF;
        output << "Stopwatch: off" << endl;
    }
    else if (!next.empty())
    {
        stopwatch_mode = StopwatchMode::NEXT;
        output << "Stopwatch: on for the next command" << endl;
    }
    else
    {
        assert(!"Impossible stopwatch mode!");
    }

    return {};
}

std::string MainProgram::print_affiliation_name(AffiliationID id, std::ostream &output, bool nl)
{
    try
    {
        if (id != NO_AFFILIATION)
        {
            auto name = ds_.get_affiliation_name(id);
            if (!name.empty())
            {
                output << name;
            }
            else
            {
                output << "*";
            }

            ostringstream retstream;
            retstream << name;
            if (nl) { output << endl; }
            return retstream.str();
        }
        else
        {
            output << "--NO_AFFILIATION--";
            if (nl) { output << endl; }
            return "";
        }
    }
    catch (NotImplemented const& e)
    {
        output << endl << "NotImplemented while printing affiliation name : " << e.what() << endl;
        std::cerr << endl << "NotImplemented while printing affiliation name : " << e.what() << endl;
        return "";
    }
}

std::string MainProgram::print_coord(Coord coord, std::ostream& output, bool nl)
{
    if (coord != NO_COORD)
    {
        output << "(" << coord.x << "," << coord.y << ")";
        ostringstream retstream;
        retstream << "(" << coord.x << "," << coord.y << ")";
        if (nl) { output << endl; }
        return retstream.str();
    }
    else
    {
        output << "(--NO_COORD--)";
        if (nl) { output << endl; }
        return "";
    }
}

string const affiliationidx = "([a-zA-Z0-9-]+)";
string const affiliationlistx = "(?:[a-zA-Z0-9-]+)";
string const publicationidx = "([0-9]+)";
string const namex = "([ a-zA-Z0-9-]+)";
string const timex = "([0-9]+)";
string const numx = "([0-9]+)";
string const optcoordx = "\\([[:space:]]*[0-9]+[[:space:]]*,[[:space:]]*[0-9]+[[:space:]]*\\)";
string const coordx = "\\([[:space:]]*([0-9]+)[[:space:]]*,[[:space:]]*([0-9]+)[[:space:]]*\\)";
string const wsx = "[[:space:]]+";


vector<MainProgram::CmdInfo> MainProgram::cmds_ =
{
    {"get_affiliation_count", "", "", &MainProgram::cmd_get_affiliation_count, &MainProgram::test_get_affiliation_count },
    {"clear_all", "", "", &MainProgram::cmd_clear_all, nullptr }, // clear all probably shouldn't be perftested since it will ... clear everything
    {"get_all_affiliations", "", "", &MainProgram::cmd_get_all_affiliations, &MainProgram::NoParListTestCmd<&Datastructures::get_all_affiliations>},
    {"add_affiliation", "AffiliationID \"Name\" (x,y)", affiliationidx+wsx+'"'+namex+'"'+wsx+coordx, &MainProgram::cmd_add_affiliation, nullptr }, // tested within each perftest, separate perftesting not necessary
    {"affiliation_info", "AffiliationID", affiliationidx, &MainProgram::cmd_affiliation_info, &MainProgram::test_affiliation_info },
    {"get_affiliations_alphabetically", "", "", &MainProgram::NoParListCmd<&Datastructures::get_affiliations_alphabetically>, &MainProgram::NoParListTestCmd<&Datastructures::get_affiliations_alphabetically> },
    {"get_affiliations_distance_increasing", "", "", &MainProgram::NoParListCmd<&Datastructures::get_affiliations_distance_increasing>,
     &MainProgram::NoParListTestCmd<&Datastructures::get_affiliations_distance_increasing> },
    {"find_affiliation_with_coord", "(x,y)", coordx, &MainProgram::cmd_find_affiliation_with_coord, &MainProgram::test_find_affiliation_with_coord },
    {"change_affiliation_coord", "AffiliationID (x,y)", affiliationidx+wsx+coordx, &MainProgram::cmd_change_affiliation_coord, &MainProgram::test_change_affiliation_coord },
    {"get_publications_after", "AffiliationID Time", affiliationidx+wsx+timex, &MainProgram::cmd_get_publications_after, &MainProgram::test_get_publications_after },
    {"add_publication", "PublicationID \"Name\" Year AffiliationID AffiliationID ...", publicationidx+wsx+'"'+namex+'"'+wsx+timex+"((?:"+wsx+affiliationlistx+")*)", &MainProgram::cmd_add_publication, nullptr }, // tested within each perftest, separate perftesting not necessary
    {"get_all_publications", "", "", &MainProgram::cmd_get_all_publications, &MainProgram::test_get_all_publications},
    {"publication_info", "PublicationID", publicationidx, &MainProgram::cmd_publication_info, &MainProgram::test_publication_info },
    {"add_reference", "PublicationID parentPublicationID", publicationidx+wsx+publicationidx, &MainProgram::cmd_add_reference, nullptr },
    {"add_affiliation_to_publication", "AffiliationID PublicationID", affiliationidx+wsx+publicationidx, &MainProgram::cmd_add_affiliation_to_publication, &MainProgram::test_add_affiliation_to_publication},
    {"get_publications", "AffiliationID", affiliationidx, &MainProgram::cmd_get_publications, &MainProgram::test_get_publications },
    {"get_all_references", "PublicationID", publicationidx, &MainProgram::cmd_get_all_references, &MainProgram::test_get_all_references },
    {"get_affiliations_closest_to", "(x,y)", coordx, &MainProgram::cmd_get_affiliations_closest_to, &MainProgram::test_affiliations_closest_to },
    {"remove_affiliation", "AffiliationID", affiliationidx, &MainProgram::cmd_remove_affiliation, &MainProgram::test_remove_affiliation },
    {"get_closest_common_parent", "PublicationID1 PublicationID2", publicationidx+wsx+publicationidx, &MainProgram::cmd_get_closest_common_parent, &MainProgram::test_get_closest_common_parent },
    {"quit", "", "", nullptr, nullptr },
    {"help", "", "", &MainProgram::help_command, nullptr },
    {"random_add", "number_of_affiliations_to_add  (minx,miny) (maxx,maxy) (coordinates optional)",
     numx+"(?:"+wsx+coordx+wsx+coordx+")?", &MainProgram::cmd_random_affiliations, &MainProgram::test_random_affiliations },
    {"read", "\"in-filename\" [silent]", "\"([-a-zA-Z0-9 ./:_]+)\"(?:"+wsx+"(silent))?", &MainProgram::cmd_read, nullptr },
    {"testread", "\"in-filename\" \"out-filename\"", "\"([-a-zA-Z0-9 ./:_]+)\""+wsx+"\"([-a-zA-Z0-9 ./:_]+)\"", &MainProgram::cmd_testread, nullptr },
    {"perftest", "cmd1[;cmd2...] timeout repeat_count n1[;n2...] (parts in [] are optional, alternatives separated by |)",
     "([0-9a-zA-Z_]+(?:;[0-9a-zA-Z_]+)*)"+wsx+numx+wsx+numx+wsx+"([0-9]+(?:;[0-9]+)*)", &MainProgram::cmd_perftest, nullptr },
    {"stopwatch", "on|off|next (alternatives separated by |)", "(?:(on)|(off)|(next))", &MainProgram::cmd_stopwatch, nullptr },
    {"random_seed", "new-random-seed-integer", numx, &MainProgram::cmd_randseed, nullptr },
    {"#", "comment text", ".*", &MainProgram::cmd_comment, nullptr },
    {"remove_publication","PublicationID",publicationidx, &MainProgram::cmd_remove_publication, &MainProgram::test_remove_publication},
    {"get_parent","PublicationID",publicationidx,&MainProgram::cmd_get_parent, &MainProgram::test_get_parent},
    {"get_referenced_by_chain","PublicationID",publicationidx,&MainProgram::cmd_get_referenced_by_chain,&MainProgram::test_get_referenced_by_chain},
    {"get_affiliations", "PublicationID", publicationidx, &MainProgram::cmd_get_affiliations, &MainProgram::test_get_affiliations},
    {"get_direct_references", "PublicationID", publicationidx, &MainProgram::cmd_get_direct_references, &MainProgram::test_get_direct_references},
    // prg2
    {"get_connected_affiliations","AffiliationID", affiliationidx, &MainProgram::cmd_get_connected_affiliations,&MainProgram::test_get_connected_affiliations},
    {"get_all_connections","","",&MainProgram::cmd_get_all_connections,&MainProgram::test_get_all_connections},
    {"get_any_path", "AffiliationID AffiliationID", affiliationidx+wsx+affiliationidx,&MainProgram::cmd_get_any_path,&MainProgram::test_get_any_path},
    // prg2 optional
    {"get_path_with_least_affiliations", "AffiliationID AffiliationID", affiliationidx+wsx+affiliationidx,&MainProgram::cmd_get_path_with_least_affiliations,&MainProgram::test_get_path_with_least_affiliations},
    {"get_path_of_least_friction", "AffiliationID AffiliationID", affiliationidx+wsx+affiliationidx,&MainProgram::cmd_get_path_of_least_friction,&MainProgram::test_get_path_of_least_friction},
    {"get_shortest_path", "AffiliationID AffiliationID", affiliationidx+wsx+affiliationidx,&MainProgram::cmd_get_shortest_path,&MainProgram::test_get_shortest_path},

};

MainProgram::CmdResult MainProgram::help_command(std::ostream& output, MatchIter /*begin*/, MatchIter /*end*/)
{
    output << "Commands:" << endl;
    for (auto& i : cmds_)
    {
        output << "  " << i.cmd << " " << i.info << endl;
    }

    return {};
}

MainProgram::CmdResult MainProgram::cmd_perftest(std::ostream& output, MatchIter begin, MatchIter end)
{
#ifdef _GLIBCXX_DEBUG
    output << "WARNING: Debug STL enabled, performance will be worse than expected (maybe also asymptotically)!" << endl;
#endif // _GLIBCXX_DEBUG

    try {
    // Note: everything below is indented too little by one indentation level! (because of try block above)

    string commandstr = *begin++;
    unsigned int timeout = convert_string_to<unsigned int>(*begin++);
    unsigned int repeat_count = convert_string_to<unsigned int>(*begin++);
    string sizes = *begin++;
    assert(begin == end && "Invalid number of parameters");

    vector<string> testcmds;
    smatch scmd;
    auto cbeg = commandstr.cbegin();
    auto cend = commandstr.cend();
    for ( ; regex_search(cbeg, cend, scmd, commands_regex_); cbeg = scmd.suffix().first)
    {
        testcmds.push_back(scmd[1]);
    }


    vector<unsigned int> init_ns;
    smatch size;
    auto sbeg = sizes.cbegin();
    auto send = sizes.cend();
    for ( ; regex_search(sbeg, send, size, sizes_regex_); sbeg = size.suffix().first)
    {
        init_ns.push_back(convert_string_to<unsigned int>(size[1]));
    }

    output << "Timeout for each N is " << timeout << " sec. " << endl;
    output << "For each N perform " << repeat_count << " random command(s) from:" << endl;

    // Initialize test functions
    vector<void(MainProgram::*)()> testfuncs;

    for (auto& i : testcmds)
    {
        auto pos = find_if(cmds_.begin(), cmds_.end(), [&i](auto const& cmd){ return cmd.cmd == i; });
        if (pos != cmds_.end() && pos->testfunc)
        {
            output << i << " ";
            testfuncs.push_back(pos->testfunc);
        }
        else
        {
            output << "(cannot test " << i << ") ";
        }
    }

    output << endl << endl;

    if (testfuncs.empty())
    {
        output << "No commands to test!" << endl;
        return {};
    }

#ifdef USE_PERF_EVENT
    output << setw(7) << "N" << " , " << setw(12) << "add (sec)" << " , " << setw(12) << "add (count)" << " , " << setw(12) << "cmds (sec)" << " , "
           << setw(12) << "cmds (count)"  << " , " << setw(12) << "total (sec)" << " , " << setw(12) << "total (count)" << endl;
#else
    output << setw(7) << "N" << " , " << setw(12) << "add (sec)" << " , " << setw(12) << "cmds (sec)" << " , "
           << setw(12) << "total (sec)" << endl;
#endif
    flush_output(output);

    auto stop = false;
    for (unsigned int n : init_ns)
    {
        if (stop) { break; }

        output << setw(7) << n << " , " << flush;

        ds_.clear_all();
        init_primes();

        Stopwatch stopwatch(true); // Use also instruction counting, if enabled
        std::unordered_set<Coord,CoordHash> exclude_list;
        std::vector<Coord> unique_coords = get_unique_coords(n,exclude_list,RANDOM_MIN_COORD,RANDOM_MAX_COORD);
        // Add random affiliations
        std::vector<Coord>::iterator start_of_range=unique_coords.begin();
        for (unsigned int i = 0; i < n / 1000; ++i,std::advance(start_of_range,1000))
        {
            std::vector<Coord> vector_slice(start_of_range,std::next(start_of_range,1000));
            stopwatch.start();
            add_random_affiliations_publications(1000,RANDOM_MIN_COORD,RANDOM_MAX_COORD,vector_slice);
            stopwatch.stop();

            if (stopwatch.elapsed() >= timeout)
            {
                output << "ADD Timeout!" << endl;
                stop = true;
                break;
            }
            if (check_stop())
            {
                output << "Stopped!" << endl;
                stop = true;
                break;
            }
        }
        if (stop) { break; }

        if (n % 1000 != 0)
        {
            std::vector<Coord> vector_slice(start_of_range,unique_coords.end());
            stopwatch.start();
            add_random_affiliations_publications(n % 1000,RANDOM_MIN_COORD,RANDOM_MAX_COORD,vector_slice);
            stopwatch.stop();
        }

#ifdef USE_PERF_EVENT
        auto addcount = stopwatch.count();
#endif
        auto addsec = stopwatch.elapsed();

#ifdef USE_PERF_EVENT
        output << setw(12) << addsec << " , " << setw(12) << addcount << " , " << flush;
#else
        output << setw(12) << addsec << " , " << flush;
#endif

        if (addsec >= timeout)
        {
            output << "ADD Timeout!" << endl;
            stop = true;
            break;
        }

        stopwatch.start();
        for (unsigned int repeat = 0; repeat < repeat_count; ++repeat)
        {
            auto cmdpos = random(testfuncs.begin(), testfuncs.end());

            (this->**cmdpos)();

            if (repeat % 10 == 0)
            {
                stopwatch.stop();
                if (stopwatch.elapsed() >= timeout)
                {
                    output << "Timeout!" << endl;
                    stop = true;
                    break;
                }
                if (check_stop())
                {
                    output << "Stopped!" << endl;
                    stop = true;
                    break;
                }
                stopwatch.start();
            }
        }
        stopwatch.stop();
        if (stop) { break; }

#ifdef USE_PERF_EVENT
        auto totalcount = stopwatch.count();
#endif
        auto totalsec = stopwatch.elapsed();

#ifdef USE_PERF_EVENT
        output << setw(12) << totalsec-addsec << " , " << setw(12) << totalcount-addcount << " , " << setw(12) << totalsec << " , " << setw(12) << totalcount;
#else
        output << setw(12) << totalsec-addsec << " , " << setw(12) << totalsec;
#endif

        output << endl;
        flush_output(output);
    }

    ds_.clear_all();
    init_primes();

    }
    catch (NotImplemented const&)
    {
        // Clean up after NotImplemented
        ds_.clear_all();
        init_primes();
        throw;
    }

#ifdef _GLIBCXX_DEBUG
    output << "WARNING: Debug STL enabled, performance will be worse than expected (maybe also asymptotically)!" << endl;
#endif // _GLIBCXX_DEBUG

    return {};
}

MainProgram::CmdResult MainProgram::cmd_comment(std::ostream& /*output*/, MatchIter /*begin*/, MatchIter /*end*/)
{
    return {};
}


bool MainProgram::command_parse_line(string inputline, ostream& output)
{

    if (inputline.empty()) { return true; }

    smatch match;
    bool matched = regex_match(inputline, match, cmds_regex_);
    if (matched)
    {
        assert(match.size() == 3);
        string cmd = match[1];
        string params = match[2];

        auto pos = find_if(cmds_.begin(), cmds_.end(), [cmd](CmdInfo const& ci) { return ci.cmd == cmd; });
        assert(pos != cmds_.end());

        smatch match2;
        bool matched2 = regex_match(params, match2, pos->param_regex);
        if (matched2)
        {
            if (pos->func)
            {
                assert(!match2.empty());

                Stopwatch stopwatch(true);
                bool use_stopwatch = (stopwatch_mode != StopwatchMode::OFF);
                // Reset stopwatch mode if only for the next command
                if (stopwatch_mode == StopwatchMode::NEXT) { stopwatch_mode = StopwatchMode::OFF; }

               TestStatus initial_status = test_status_;
               test_status_ = TestStatus::NOT_RUN;

                if (use_stopwatch)
                {
                    stopwatch.start();
                }

                CmdResult result;
                try
                {
                    result = (this->*(pos->func))(output, ++(match2.begin()), match2.end());
                }
                catch (NotImplemented const& e)
                {
                    output << endl << "NotImplemented from cmd " << pos->cmd << " : " << e.what() << endl;
                    std::cerr << endl << "NotImplemented from cmd " << pos->cmd << " : " << e.what() << endl;
                }

                if (use_stopwatch)
                {
                    stopwatch.stop();
                }

                switch (result.first)
                {
                    case ResultType::NOTHING:
                    {
                        break;
                    }
                    case ResultType::IDLIST:
                    {
                        auto& [publications, affiliations] = std::get<CmdResultIDs>(result.second);
                        if (affiliations.size() == 1 && affiliations.front() == NO_AFFILIATION)
                        {
                            output << "Failed (NO_AFFILIATION returned)!" << std::endl;
                        }
                        else
                        {
                            if (!affiliations.empty())
                            {
                                if (affiliations.size() == 1) { output << "Affiliation:" << std::endl; }
                                else { output << "Affiliations:" << std::endl; }

                                unsigned int num = 0;
                                for (AffiliationID& id : affiliations)
                                {
                                    ++num;
                                    if (affiliations.size() > 1) { output << num << ". "; }
                                    else { output << "   "; }
                                    print_affiliation(id, output);
                                }
                            }
                        }

                        if (publications.size() == 1 && publications.front() == NO_PUBLICATION)
                        {
                            output << "Failed (NO_PUBLICATION returned)!" << std::endl;
                        }
                        else
                        {
                            if (!publications.empty())
                            {
                                if (publications.size() == 1) { output << "Publication:" << std::endl; }
                                else { output << "Publications:" << std::endl; }

                                unsigned int num = 0;
                                for (PublicationID id : publications)
                                {
                                    ++num;
                                    if (publications.size() > 1) { output << num << ". "; }
                                    else { output << "   "; }
                                    print_publication(id, output);
                                }
                            }
                        }
                        break;
                    }
                    case ResultType::ROUTE:
                    {
                        auto& route = std::get<CmdResultRoute>(result.second);
                        if (!route.empty())
                        {
                            if (route.size() == 1 && get<0>(route.front()) == NO_AFFILIATION)
                            {
                                output << "Failed (...NO_AFFILIATION... returned)!" << std::endl;
                            }
                            else
                            {
                                unsigned int num = 1;
                                for (auto& r : route)
                                {
                                    auto [affiliationid1, weight, affiliationid2, dist] = r;
                                    output << num << ". ";
                                    if (affiliationid1 != NO_AFFILIATION)
                                    {
                                        print_affiliation_brief(affiliationid1, output, false);
                                    }
                                    if (affiliationid2 != NO_AFFILIATION)
                                    {
                                        output << " -> ";
                                        print_affiliation_brief(affiliationid2, output, false);
                                    }
                                    if (weight != NO_WEIGHT)
                                    {
                                        output << " (weighted " << weight << ")";
                                    }
                                    if (dist != NO_DISTANCE)
                                    {
                                        output << " (distance " << dist << ")";
                                    }
                                    output << endl;

                                    ++num;
                                }
                            }
                        }
                        break;
                    }
                    case ResultType::CONNECTIONLIST:{
                        auto& list = std::get<ConnectionList>(result.second);
                        unsigned int num = 1;

                        std::for_each(list.begin(),list.end(),[&output,&num,this](auto& connection){
                            output << num++ << ". ";
                            print_affiliation_brief(connection.aff1 ,output,false);
                            output << " -> ";
                            print_affiliation_brief(connection.aff2, output, false);
                            output <<" (weighted "<<connection.weight<<")"<<endl;
                        });
                        break;
                    }
                    case ResultType::NEIGHBOURLIST:{
                        auto& list = std::get<ConnectionList>(result.second);
                        auto source_id = list.front().aff1;
                        output << "All connected affiliations from ";
                        print_affiliation_brief(source_id,output,false);
                        output << endl;
                        unsigned int num = 1;
                        std::for_each(list.begin(),list.end(),[&output,&num,this](auto& connection){
                            output << num++ << ". ";
                            print_affiliation_brief(connection.aff2, output, false);
                            output <<" (weighted "<<connection.weight<<")"<<endl;
                        });
                        break;
                    }
                    default:
                    {
                        assert(false && "Unsupported result type!");
                    }
                }

                if (result != prev_result)
                {
                    prev_result = move(result);
                    view_dirty = true;
                }

                if (use_stopwatch)
                {
                    output << "Command '" << cmd << "': " << stopwatch.elapsed() << " sec";
#ifdef USE_PERF_EVENT
                    auto totalcount = stopwatch.count();
                    output << ", cmds (count): " << totalcount;
#endif
                    output << endl;
                }

                if (test_status_ != TestStatus::NOT_RUN)
                {
                    output << "Testread-tests have been run, " << ((test_status_ == TestStatus::DIFFS_FOUND) ? "differences found!" : "no differences found.") << endl;
                }
                if (test_status_ == TestStatus::NOT_RUN || (test_status_ == TestStatus::NO_DIFFS && initial_status == TestStatus::DIFFS_FOUND))
                {
                    test_status_ = initial_status;
                }
            }
            else
            { // No function to run = quit command
                return false;
            }
        }
        else
        {
            output << "Invalid parameters for command '" << cmd << "'!" << endl;
        }
    }
    else
    {
        output << "Unknown command!" << endl;
    }

    return true; // Signal continuing
}

void MainProgram::command_parser(istream& input, ostream& output, PromptStyle promptstyle)
{
    string line;
    do
    {
        output << PROMPT;
        getline(input, line, '\n');

        if (promptstyle != PromptStyle::NO_ECHO)
        {
            output << line << endl;
        }

        if (!input) { break; }

        bool cont = command_parse_line(line, output);
        view_dirty = false; // No need to keep track of individual result changes
        if (!cont) { break; }
    }
    while (input);

    view_dirty = true; // To be safe, assume that results have been changed
}

void MainProgram::setui(MainWindow* ui)
{
    ui_ = ui;
}

#ifdef GRAPHICAL_GUI
void MainProgram::flush_output(std::ostream& output)
{
    if (ui_)
    {
        if (auto soutput = dynamic_cast<ostringstream*>(&output))
        {
            ui_->output_text(*soutput);
        }
    }
}
#else
void MainProgram::flush_output(std::ostream& /*output*/)
{
}
#endif

bool MainProgram::check_stop() const
{
#ifdef GRAPHICAL_GUI
    if (ui_)
    {
        return ui_->check_stop_pressed();
    }
#endif
    return false;
}

std::array<unsigned long int, 20> const MainProgram::primes1{4943,   4951,   4957,   4967,   4969,   4973,   4987,   4993,   4999,   5003,
                                                             5009,   5011,   5021,   5023,   5039,   5051,   5059,   5077,   5081,   5087};
std::array<unsigned long int, 20> const MainProgram::primes2{81031,  81041,  81043,  81047,  81049,  81071,  81077,  81083,  81097,  81101,
                                                             81119,  81131,  81157,  81163,  81173,  81181,  81197,  81199,  81203,  81223};

MainProgram::MainProgram()
{
    rand_engine_.seed(time(nullptr));

    init_primes();
    init_regexs();
}

int MainProgram::mainprogram(int argc, char* argv[])
{
    vector<string> args(argv, argv+argc);

    if (args.size() < 1 || args.size() > 2)
    {
        cerr << "Usage: " + ((args.size() > 0) ? args[0] : "<program name>") + " [<command file>]" << endl;
        return EXIT_FAILURE;
    }

    MainProgram mainprg;

    if (args.size() == 2 && args[1] != "--console")
    {
        string filename = args[1];
        ifstream input(filename);
        if (input)
        {
            mainprg.command_parser(input, cout, MainProgram::PromptStyle::NORMAL);
        }
        else
        {
            cout << "Cannot open file '" << filename << "'!" << endl;
        }
    }
    else
    {
        mainprg.command_parser(cin, cout, MainProgram::PromptStyle::NO_ECHO);
    }

    cerr << "Program ended normally." << endl;
    if (mainprg.test_status_ == TestStatus::DIFFS_FOUND)
    {
        return EXIT_FAILURE;
    }
    else
    {
        return EXIT_SUCCESS;
    }
}

void MainProgram::init_primes()
{
    // Initialize id generator
    prime1_ = primes1[random<int>(0, primes1.size())];
    prime2_ = primes2[random<int>(0, primes2.size())];
    random_affiliations_added_ = 0;
    random_publications_added_ = 0;
}

Name MainProgram::n_to_name(unsigned long n)
{
    unsigned long int hash = prime1_*n + prime2_;
    string name;

    while (hash > 0)
    {
        auto hexnum = hash % 26;
        hash /= 26;
        name.push_back('a'+hexnum);
    }

    return name;
}

AffiliationID MainProgram::n_to_affiliationid(unsigned long n)
{
 std::ostringstream ostr;
 ostr << "A" << n;
 return ostr.str();
}

PublicationID MainProgram::n_to_publicationid(unsigned long n)
{
    return n;
}

Coord MainProgram::n_to_coord(unsigned long n)
{
    unsigned long int hash = prime1_ * n + prime2_;
    hash = hash ^ (hash + 0x9e3779b9 + (hash << 6) + (hash >> 2)); // :-P

    return {static_cast<int>(hash % 1000), static_cast<int>((hash/1000) % 1000)};
}

void MainProgram::init_regexs()
{
    // Create regex <whitespace>(cmd1|cmd2|...)<whitespace>(.*)
    string cmds_regex_str = "[[:space:]]*(";
    bool first = true;
    for (auto& cmd : cmds_)
    {
        cmds_regex_str += (first ? "" : "|") + cmd.cmd;
        first = false;

        cmd.param_regex = regex(cmd.param_regex_str+"[[:space:]]*", std::regex_constants::ECMAScript | std::regex_constants::optimize);
    }
    cmds_regex_str += ")(?:[[:space:]]*$|"+wsx+"(.*))";
    cmds_regex_ = regex(cmds_regex_str, std::regex_constants::ECMAScript | std::regex_constants::optimize);
    coords_regex_ = regex(coordx+"[[:space:]]?", std::regex_constants::ECMAScript | std::regex_constants::optimize);
    affil_regex_ = regex(affiliationidx+"[[:space:]]?", std::regex_constants::ECMAScript | std::regex_constants::optimize); // TODO test this one more intensively
    times_regex_ = regex(wsx+"([0-9][0-9]):([0-9][0-9]):([0-9][0-9])", std::regex_constants::ECMAScript | std::regex_constants::optimize);
    commands_regex_ = regex("([0-9a-zA-Z_]+);?", std::regex_constants::ECMAScript | std::regex_constants::optimize);
    sizes_regex_ = regex(numx+";?", std::regex_constants::ECMAScript | std::regex_constants::optimize);
}
