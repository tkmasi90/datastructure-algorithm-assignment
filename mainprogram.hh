// The main program (provided by the course), COMP.CS.300
//
// DO ****NOT**** EDIT THIS FILE!
// (Preferably do not edit this even temporarily. And if you still decide to do so
//  (for debugging, for example), DO NOT commit any changes to git, but revert all
//  changes later. Otherwise you won't be able to get any updates/fixes to this
//  file from git!)
#ifndef MAINPROGRAM_HH
#define MAINPROGRAM_HH


#ifdef QT_CORE_LIB
#ifndef NO_GRAPHICAL_GUI
#define GRAPHICAL_GUI
#endif
#endif


#include <string>
#include <random>
#include <regex>
#include <chrono>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <array>
#include <functional>
#include <utility>
#include <variant>
#include <bitset>
#include <cassert>
#include <cstring>
#include <unordered_set>

#include "datastructures.hh"

// default max and min values for perftesting and random add, may be subject to change

const Coord RANDOM_MIN_COORD = {0,0};
const Coord RANDOM_MAX_COORD = {10000,10000};
const Year RANDOM_MIN_YEAR = 0;
const Year RANDOM_MAX_YEAR = 9998;

const double ROOT_BIAS_MULTIPLIER = 0.05;
const double LEAF_BIAS_MULTIPLIER = 0.5;

class MainWindow; // In case there's UI

class MainProgram
{
public:
    MainProgram();


    class Stopwatch;

    enum class PromptStyle { NORMAL, NO_ECHO, NO_NESTING };
    enum class TestStatus { NOT_RUN, NO_DIFFS, DIFFS_FOUND };

    bool command_parse_line(std::string input, std::ostream& output);
    void command_parser(std::istream& input, std::ostream& output, PromptStyle promptstyle);

    void setui(MainWindow* ui);

    void flush_output(std::ostream& output);
    bool check_stop() const;

    static int mainprogram(int argc, char* argv[]);

private:
    Datastructures ds_;
    MainWindow* ui_ = nullptr;

    static std::string const PROMPT;

    std::minstd_rand rand_engine_;

    static std::array<unsigned long int, 20> const primes1;
    static std::array<unsigned long int, 20> const primes2;
    unsigned long int prime1_ = 0; // Will be initialized to random value from above
    unsigned long int prime2_ = 0; // Will be initialized to random value from above
    unsigned long int random_affiliations_added_ = 0; // Counter for random affiliations added
    unsigned long int random_publications_added_ = 0; // Counter for random publications added
    void init_primes();
    Name n_to_name(unsigned long int n);
    AffiliationID n_to_affiliationid(unsigned long int n);
    PublicationID n_to_publicationid(unsigned long int n);
    Coord n_to_coord(unsigned long int n);


    enum class StopwatchMode { OFF, ON, NEXT };
    StopwatchMode stopwatch_mode = StopwatchMode::OFF;

    enum class ResultType { NOTHING, IDLIST, ROUTE, CONNECTIONLIST, NEIGHBOURLIST};
    using CmdResultIDs = std::pair<std::vector<PublicationID>, std::vector<AffiliationID>>;
    using CmdResultRoute = std::vector<std::tuple<AffiliationID, Weight, AffiliationID, Distance>>;
    using ConnectionList = std::vector<Connection>;
    using CmdResult = std::pair<ResultType, std::variant<CmdResultIDs, CmdResultRoute, ConnectionList>>;
    CmdResult prev_result;
    bool view_dirty = true;

    TestStatus test_status_ = TestStatus::NOT_RUN;

    using MatchIter = std::smatch::const_iterator;
    struct CmdInfo
    {
        std::string cmd;
        std::string info;
        std::string param_regex_str;
        CmdResult(MainProgram::*func)(std::ostream& output, MatchIter begin, MatchIter end);
        void(MainProgram::*testfunc)();
        std::regex param_regex = {};
    };
    static std::vector<CmdInfo> cmds_;
    // Regex objects and their initialization
    std::regex cmds_regex_;
    std::regex coords_regex_;
    std::regex affil_regex_;
    std::regex times_regex_;
    std::regex commands_regex_;
    std::regex sizes_regex_;
    void init_regexs();


    CmdResult cmd_get_affiliation_count(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_clear_all(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_get_all_affiliations(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_add_affiliation(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_affiliation_info(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_find_affiliation_with_coord(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_change_affiliation_coord(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_get_publications_after(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_add_publication(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_get_all_publications(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_publication_info(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_add_reference(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_add_affiliation_to_publication(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_get_publications(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_get_all_references(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_get_affiliations_closest_to(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_remove_affiliation(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_get_closest_common_parent(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_remove_publication(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_get_parent(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_get_referenced_by_chain(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_get_direct_references(std::ostream& output, MatchIter begin, MatchIter end);

    CmdResult help_command(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_randseed(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_random_affiliations(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_read(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_testread(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_stopwatch(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_perftest(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_comment(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_get_affiliations(std::ostream& output, MatchIter begin, MatchIter end);
    // PRG2 command functions
    CmdResult cmd_get_connected_affiliations(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_get_all_connections(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_get_any_path(std::ostream& output, MatchIter begin, MatchIter end);
    // PRG2 optional
    CmdResult cmd_get_path_with_least_affiliations(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_get_path_of_least_friction(std::ostream& output, MatchIter begin, MatchIter end);
    CmdResult cmd_get_shortest_path(std::ostream& output, MatchIter begin, MatchIter end);

    // random ids for perftest
    AffiliationID random_affiliation();
    PublicationID random_publication();

    // biased random ids for some perftest
    PublicationID random_root_publication();
    PublicationID random_leaf_publication();

    void test_get_functions(AffiliationID id);
    void test_affiliation_info();
    void test_find_affiliation_with_coord();
    void test_change_affiliation_coord();
    void test_get_publications_after();
    void test_publication_info();
    void test_get_publications();
    void test_get_all_references();
    void test_affiliations_closest_to();
    void test_remove_affiliation();
    void test_get_closest_common_parent();
    void test_random_affiliations();
    void test_remove_publication();
    void test_get_parent();
    void test_get_referenced_by_chain();
    void test_get_direct_references();
    void test_get_affiliations();
    void test_get_affiliation_count();
    void test_get_all_publications();
    void test_add_affiliation_to_publication();
    // prg2
    void test_get_connected_affiliations();
    void test_get_all_connections();
    void test_get_any_path();
    // prg2 optional
    void test_get_path_with_least_affiliations();
    void test_get_path_of_least_friction();
    void test_get_shortest_path();


    inline Coord get_random_coords(const Coord min = RANDOM_MIN_COORD, const Coord max = RANDOM_MAX_COORD);
    inline Year get_random_year(const Year min = RANDOM_MIN_YEAR, const Year max = RANDOM_MAX_YEAR);
    std::vector<Coord> get_unique_coords(const unsigned int n,const std::unordered_set<Coord,CoordHash>& exclude_list,const Coord min=RANDOM_MIN_COORD,const Coord max=RANDOM_MAX_COORD);
    void add_random_affiliations_publications(unsigned int size, Coord min = RANDOM_MIN_COORD, Coord max = RANDOM_MAX_COORD,const std::vector<Coord>& coordinates={});
    Distance calc_distance(Coord c1, Coord c2);
    std::string print_affiliation(AffiliationID id, std::ostream& output, bool nl = true);
    std::string print_affiliation_brief(AffiliationID id, std::ostream& output, bool nl = true);
    std::string print_publication(PublicationID id, std::ostream& output, bool nl = true);
    std::string print_affiliation_name(AffiliationID id, std::ostream& output, bool nl = true);
    std::string print_coord(Coord coord, std::ostream& output, bool nl = true);

    template <typename Type>
    Type random(Type start, Type end);
    template <typename To>
    static To convert_string_to(std::string from);
    template <typename From>
    static std::string convert_to_string(From from);

    template<AffiliationID(Datastructures::*MFUNC)()>
    CmdResult NoParAffiliationCmd(std::ostream& output, MatchIter begin, MatchIter end);

    template<std::vector<AffiliationID>(Datastructures::*MFUNC)()>
    CmdResult NoParListCmd(std::ostream& output, MatchIter begin, MatchIter end);

    template<AffiliationID(Datastructures::*MFUNC)()>
    void NoParAffiliationTestCmd();

    template<std::vector<AffiliationID>(Datastructures::*MFUNC)()>
    void NoParListTestCmd();

    friend class MainWindow;
};

template <typename Type>
Type MainProgram::random(Type start, Type end)
{
    auto range = end-start;
    assert(range != 0 && "random() with zero range!");

    auto num = std::uniform_int_distribution<unsigned long int>(0, range-1)(rand_engine_);

    return static_cast<Type>(start+num);
}

template <typename To>
To MainProgram::convert_string_to(std::string from)
{
    std::istringstream istr(from);
    To result;
    istr >> std::noskipws >> result;
    if (istr.fail() || !istr.eof())
    {
        throw std::invalid_argument("Cannot convert string to required type");
    }
    return result;
}

template <typename From>
std::string MainProgram::convert_to_string(From from)
{
    std::ostringstream ostr;
    ostr << from;
    if (ostr.fail())
    {
        throw std::invalid_argument("Cannot convert type to string");
    }
    return ostr.str();
}

template<AffiliationID(Datastructures::*MFUNC)()>
MainProgram::CmdResult MainProgram::NoParAffiliationCmd(std::ostream& /*output*/, MatchIter /*begin*/, MatchIter /*end*/)
{
    auto result = (ds_.*MFUNC)();
    return {ResultType::IDLIST, CmdResultIDs{{}, {result}}};
}

template<std::vector<AffiliationID>(Datastructures::*MFUNC)()>
MainProgram::CmdResult MainProgram::NoParListCmd(std::ostream& /*output*/, MatchIter /*begin*/, MatchIter /*end*/)
{
    auto result = (ds_.*MFUNC)();
    return {ResultType::IDLIST, CmdResultIDs{{}, result}};
}

template<AffiliationID(Datastructures::*MFUNC)()>
void MainProgram::NoParAffiliationTestCmd()
{
    (ds_.*MFUNC)();
}

template<std::vector<AffiliationID>(Datastructures::*MFUNC)()>
void MainProgram::NoParListTestCmd()
{
    (ds_.*MFUNC)();
}


#ifdef USE_PERF_EVENT
extern "C"
{
#include <unistd.h>
#include <sys/ioctl.h>
#include <linux/perf_event.h>
#include <asm/unistd.h>

static long
perf_event_open(struct perf_event_attr *hw_event, pid_t pid,
                int cpu, int group_fd, unsigned long flags)
{
    int ret;

    ret = syscall(__NR_perf_event_open, hw_event, pid, cpu,
                   group_fd, flags);
    return ret;
}
}
#endif

class MainProgram::Stopwatch
{
public:
    using Clock = std::chrono::high_resolution_clock;

    Stopwatch(bool use_counter = false) : use_counter_(use_counter)
    {
#ifdef USE_PERF_EVENT
        if (use_counter_)
        {
            memset(&pe_, 0, sizeof(pe_));
            pe_.type = PERF_TYPE_HARDWARE;
            pe_.size = sizeof(pe_);
            pe_.config = PERF_COUNT_HW_INSTRUCTIONS;
            pe_.disabled = 1;
            pe_.exclude_kernel = 1;
            pe_.exclude_hv = 1;

            fd_ = perf_event_open(&pe_, 0, -1, -1, 0);
            if (fd_ == -1) {
                throw "Couldn't open perf events!";
            }
        }
#endif
        reset();
    }

    ~Stopwatch()
    {
#ifdef USE_PERF_EVENT
        if (use_counter_)
        {
            close(fd_);
        }
#endif
    }

    void start()
    {
        running_ = true;
        starttime_ = Clock::now();
#ifdef USE_PERF_EVENT
        if (use_counter_)
        {
            ioctl(fd_, PERF_EVENT_IOC_RESET, 0);
            read(fd_, &startcount_, sizeof(startcount_));
            ioctl(fd_, PERF_EVENT_IOC_ENABLE, 0);
        }
#endif
    }

    void stop()
    {
        running_ = false;
#ifdef USE_PERF_EVENT
        if (use_counter_)
        {
            ioctl(fd_, PERF_EVENT_IOC_DISABLE, 0);
            long long c; read(fd_, &c, sizeof(c));
            counter_ += (c - startcount_);
        }
#endif
        elapsed_ += (Clock::now() - starttime_);
    }

    void reset()
    {
        running_ = false;
#ifdef USE_PERF_EVENT
        if (use_counter_)
        {
            ioctl(fd_, PERF_EVENT_IOC_DISABLE, 0);
            ioctl(fd_, PERF_EVENT_IOC_RESET, 0);
            counter_ = 0;
        }
#endif
        elapsed_ = elapsed_.zero();
    }

    double elapsed()
    {
        if (!running_)
        {
            return static_cast<std::chrono::duration<double>>(elapsed_).count();
        }
        else
        {
            auto total = elapsed_ + (Clock::now() - starttime_);
            return static_cast<std::chrono::duration<double>>(total).count();
        }
    }

#ifdef USE_PERF_EVENT
    long long count()
    {
        if (use_counter_)
        {
            if (!running_)
            {
                return counter_;
            }
            else
            {
                long long c;
                read(fd_, &c, sizeof(c));
                return counter_ + (c - startcount_);
            }
        }
        else
        {
            assert(!"perf_event not enabled during StopWatch creation!");
        }
    }
#endif

private:
    std::chrono::time_point<Clock> starttime_;
    Clock::duration elapsed_ = Clock::duration::zero();
    bool running_ = false;

    bool use_counter_;
#ifdef USE_PERF_EVENT
    struct perf_event_attr pe_;
    long long startcount_ = 0;
    long long counter_ = 0;
    int fd_ = 0;
#endif
};


#endif // MAINPROGRAM_HH
