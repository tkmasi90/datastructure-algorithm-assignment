// The main program (provided by the course), COMP.CS.300
//
// DO ****NOT**** EDIT THIS FILE!
// (Preferably do not edit this even temporarily. And if you still decide to do so
//  (for debugging, for example), DO NOT commit any changes to git, but revert all
//  changes later. Otherwise you won't be able to get any updates/fixes to this
//  file from git!)

#include "mainprogram.hh"

// Qt generated main window code

#ifdef GRAPHICAL_GUI
#include <QCoreApplication>
#include <QFileDialog>
#include <QDir>
#include <QFont>
#include <QGraphicsItem>
#include <QColor>
#include <QBrush>
#include <QPen>
#include <QGraphicsItem>
#include <QVariant>

#include <string>
using std::string;

#include <sstream>
using std::ostringstream;

#include <iostream>

#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <utility>
#include <tuple>

#include <cassert>

#include "mainwindow.hh"
#include "ui_mainwindow.h"
#if __has_include("worldmap/worldmap.hh")
#include "worldmap/worldmap.hh"
#endif
// Needed to be able to store AffiliationID in QVariant (in QGraphicsItem)
Q_DECLARE_METATYPE(AffiliationID)
// The same for PublicationIDs
Q_DECLARE_METATYPE(PublicationID)
// The same for Coords (currently a pair of ints)
Q_DECLARE_METATYPE(Coord)

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    mainprg_.setui(this);

    // Execute line
    connect(ui->execute_button, &QPushButton::pressed, this, &MainWindow::execute_line);

    // Line edit
    connect(ui->lineEdit, &QLineEdit::returnPressed, this, &MainWindow::execute_line);

    // File selection
    connect(ui->file_button, &QPushButton::pressed, this, &MainWindow::select_file);

    // Command selection
    // !!!!! Sort commands in alphabetical order (should not be done here, but is)
    std::sort(mainprg_.cmds_.begin(), mainprg_.cmds_.end(), [](auto const& l, auto const& r){ return l.cmd < r.cmd; });
    for (auto& cmd : mainprg_.cmds_)
    {
        ui->cmd_select->addItem(QString::fromStdString(cmd.cmd));
    }
    connect(ui->cmd_select, static_cast<void(QComboBox::*)(int)>(&QComboBox::activated), this, &MainWindow::cmd_selected);

    // Number selection
    for (auto i = 0; i <= 20; ++i)
    {
        ui->number_select->addItem(QString("%1").arg(i));
    }
    connect(ui->number_select, static_cast<void(QComboBox::*)(QString const&)>(&QComboBox::currentTextChanged), this, &MainWindow::number_selected);

    // Output box
    QFont monofont = QFontDatabase::systemFont(QFontDatabase::FixedFont);
    monofont.setStyleHint(QFont::TypeWriter);
    ui->output->setFont(monofont);

    // Initialize graphics scene & view
    gscene_ = new QGraphicsScene(this);
    ui->graphics_view->setScene(gscene_);
    ui->graphics_view->resetTransform();
    ui->graphics_view->setBackgroundBrush(Qt::black);

    // Selecting graphics items by mouse
    connect(gscene_, &QGraphicsScene::selectionChanged, this, &MainWindow::scene_selection_change);

    // Zoom slider changes graphics view scale
    connect(ui->zoom_plus, &QToolButton::clicked, this, [this]{ this->ui->graphics_view->scale(1.1, 1.1); });
    connect(ui->zoom_minus, &QToolButton::clicked, this, [this]{ this->ui->graphics_view->scale(1/1.1, 1/1.1); });
    connect(ui->zoom_1, &QToolButton::clicked, this, [this]{ this->ui->graphics_view->resetTransform(); });
    connect(ui->zoom_fit, &QToolButton::clicked, this, &MainWindow::fit_view);

    // Changing checkboxes updates view
    connect(ui->affiliations_checkbox, &QCheckBox::clicked, this, &MainWindow::update_view);
    connect(ui->affiliationnames_checkbox, &QCheckBox::clicked, this, &MainWindow::update_view);
    connect(ui->publications_checkbox, &QCheckBox::clicked, this, &MainWindow::update_view);
    connect(ui->connections_checkbox, &QCheckBox::clicked, this, &MainWindow::update_view);
    connect(ui->world_map_checkbox, &QCheckBox::clicked, this, &MainWindow::update_view);

    // Unchecking affiliations checkbox disables affiliation names checkbox
    connect(ui->affiliations_checkbox, &QCheckBox::clicked, this,
            [this]{ this->ui->affiliationnames_checkbox->setEnabled(this->ui->affiliations_checkbox->isChecked()); });

    // Unchecking the connections checkbox disables the min weight spinbox
    connect(ui->connections_checkbox, &QCheckBox::clicked, this, [this]{
        const bool spinbox_enabled = ui->connections_checkbox->isChecked();
        this->ui->min_weight_spinbox->setEnabled(spinbox_enabled);
        this->ui->min_weight_label->setEnabled(spinbox_enabled);
    });
    connect(ui->min_weight_spinbox,static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged),this,&MainWindow::update_view);

    // Changing font or points scale updates view
    connect(ui->fontscale, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, &MainWindow::update_view);
    connect(ui->pointscale, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, &MainWindow::update_view);

    // Clear input button
    connect(ui->clear_input_button, &QPushButton::clicked, this, &MainWindow::clear_input_line);

    // Stop button
    connect(ui->stop_button, &QPushButton::clicked, this, [this](){ this->stop_pressed_ = true; });

    // disable the worldmap checkbox if no worldmap data exists
    #ifndef WORLDMAP_HH
    ui->world_map_checkbox->setDisabled(true);
    #endif

    clear_input_line();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::update_view()
{
    std::unordered_set<std::string> errorset;
    try
    {
        gscene_->clear();
        auto pointscale = ui->pointscale->value();
        auto fontscale = ui->fontscale->value();

        // Collect the result of previous operation
        std::unordered_map<AffiliationID, std::string> result_affiliations;
        std::unordered_map<PublicationID, std::string> result_publications;
        std::vector<Connection> result_connections;
        std::map<std::pair<AffiliationID, AffiliationID>, string> result_routes; // Should be unordered_map, but pair lacks hash...
        switch (mainprg_.prev_result.first)
        {
        case MainProgram::ResultType::IDLIST:
            {
                // Copy the id vectors to the result maps
                auto& prev_result = std::get<MainProgram::CmdResultIDs>(mainprg_.prev_result.second);
                int i = 0;
                std::for_each(prev_result.second.begin(), prev_result.second.end(),
                              [&result_affiliations, &i](auto id){ result_affiliations[id] += MainProgram::convert_to_string(++i)+". "; });
                i = 0;
                std::for_each(prev_result.first.begin(), prev_result.first.end(),
                              [&result_publications, &i](auto id){ result_publications[id] += MainProgram::convert_to_string(++i)+". "; });
            }
            break;
        case MainProgram::ResultType::CONNECTIONLIST:{
                auto prev_result = std::get<MainProgram::ConnectionList>(mainprg_.prev_result.second);
                if (!prev_result.empty()){
                    result_connections = std::move(prev_result);
                }
            }
            break;
        case MainProgram::ResultType::NEIGHBOURLIST:{
                auto prev_result = std::get<MainProgram::ConnectionList>(mainprg_.prev_result.second);
                if (!prev_result.empty()){
                    result_connections = std::move(prev_result);
                }
            }
            break;
        case MainProgram::ResultType::ROUTE:{
                auto prev_result = std::get<MainProgram::CmdResultRoute>(mainprg_.prev_result.second);
                result_connections.reserve(prev_result.size());
                if (!prev_result.empty()){
                    std::transform(prev_result.begin(),prev_result.end(),std::back_inserter(result_connections),[](auto& tuple){
                        return Connection{std::get<0>(tuple),std::get<2>(tuple),std::get<1>(tuple)};
                    });
                }
            }
            break;
        case MainProgram::ResultType::NOTHING:
            break;
        default:
            assert(!"Unhandled result type in update_view()!");
        }

        auto affiliations = mainprg_.ds_.get_all_affiliations();
        if (affiliations.size() == 1 && affiliations.front() == NO_AFFILIATION)
        {
            errorset.insert("get_all_affiliations() returned error {NO_AFFILIATION}");
        }

        if (ui->affiliations_checkbox->isChecked())
        {
            auto affiliations = mainprg_.ds_.get_all_affiliations();
            if (affiliations.size() == 1 && affiliations.front() == NO_AFFILIATION)
            {
                errorset.insert("Error from GUI: get_all_affiliations() returned error {NO_AFFILIATION}");
                affiliations.clear(); // Clear the affiliations so that no more errors are caused by NO_AFFILIATION
            }

            for (auto& affiliationid : affiliations)
            {
                QColor affiliationcolor = Qt::gray;
                QColor namecolor = Qt::cyan;
                QColor affiliationborder = Qt::gray;
                int affiliationzvalue = 1;

                try
                {
                    if (affiliationid != NO_AFFILIATION)
                    {
                        auto xy = mainprg_.ds_.get_affiliation_coord(affiliationid);
                        auto [x,y] = xy;
                        if (x == NO_VALUE || y == NO_VALUE)
                        {
                            errorset.insert("get_affiliation_coordinates() returned error NO_COORD/NO_VALUE");
                        }

                        if (x == NO_VALUE || y == NO_VALUE)
                        {
                            x = 0; y = 0;
                            affiliationcolor = Qt::magenta;
                            namecolor = Qt::magenta;
                            affiliationzvalue = 30;
                        }

                        string prefix;
                        auto res_place = result_affiliations.find(affiliationid);
                        if (res_place != result_affiliations.end())
                        {
                            if (result_affiliations.size() > 1) { prefix = res_place->second; }
                            namecolor = Qt::red;
                            affiliationborder = Qt::red;
                            affiliationzvalue = 2;
                        }

                        if (ui->affiliations_checkbox->isChecked())
                        {
                            auto groupitem = gscene_->createItemGroup({});
                            groupitem->setFlag(QGraphicsItem::ItemIsSelectable);
                            groupitem->setData(0, QVariant::fromValue(affiliationid));

                            QPen placepen(affiliationborder);
                            placepen.setWidth(0); // Cosmetic pen
                            double publication_scale = std::max(1.0,1.0+std::log10(mainprg_.ds_.get_publications(affiliationid).size()));
                            auto dotitem = gscene_->addEllipse(-4*pointscale*publication_scale, -4*pointscale*publication_scale, 8*pointscale*publication_scale, 8*pointscale*publication_scale,
                                                               placepen, QBrush(affiliationcolor));
                            dotitem->setFlag(QGraphicsItem::ItemIgnoresTransformations);
                            groupitem->addToGroup(dotitem);

                            // Draw place names
                            string label = prefix;
                            if (ui->affiliationnames_checkbox->isChecked())
                            {
                                try
                                {
                                    auto name = mainprg_.ds_.get_affiliation_name(affiliationid);
                                    if (name == NO_NAME)
                                    {
                                        errorset.insert("get_affiliation_name() returned error NO_NAME");
                                    }

                                    label += name;
                                }
                                catch (NotImplemented const& e)
                                {
                                    errorset.insert(std::string("NotImplemented while updating graphics: ") + e.what());
                                    std::cerr << std::endl << "NotImplemented while updating graphics: " << e.what() << std::endl;
                                }
                            }

                            if (!label.empty())
                            {
                                // Create extra item group to be able to set ItemIgnoresTransformations on the correct level (addSimpleText does not allow
                                // setting initial coordinates in item coordinates
                                auto textgroupitem = gscene_->createItemGroup({});
                                auto textitem = gscene_->addSimpleText(QString::fromStdString(label));
                                auto font = textitem->font();
                                font.setPointSizeF(font.pointSizeF()*fontscale);
                                textitem->setFont(font);
                                textitem->setBrush(QBrush(namecolor));
                                textitem->setPos(-textitem->boundingRect().width()/2, -4*pointscale - textitem->boundingRect().height());
                                textgroupitem->addToGroup(textitem);
                                textgroupitem->setFlag(QGraphicsItem::ItemIgnoresTransformations);
                                groupitem->addToGroup(textgroupitem);
                            }

                            groupitem->setPos(20*x, -20*y);
                            groupitem->setZValue(affiliationzvalue);
                        }
                    }
                }
                catch (NotImplemented const& e)
                {
                    errorset.insert(std::string("NotImplemented while updating graphics: ") + e.what());
                    std::cerr << std::endl << "NotImplemented while updating graphics: " << e.what() << std::endl;
                }
            }
        }
        if (ui->publications_checkbox->isChecked())
        {
            try
            {
                auto publicationids = mainprg_.ds_.all_publications();
                if (publicationids.size() == 1 && publicationids.front() == NO_PUBLICATION)
                {
                    errorset.insert("all_publications() returned error {NO_PUBLICATION}");
                }

                for (auto publicationid : publicationids)
                {
                    QColor publicationcolor = Qt::blue;
                    int publicationzvalue = -3;

                    try
                    {
                        if (publicationid != NO_PUBLICATION)
                        {
                            if (result_publications.find(publicationid) != result_publications.end())
                            {
                                publicationcolor = Qt::green;
                                publicationzvalue = -2;
                            }
                            auto affiliations = mainprg_.ds_.get_affiliations(publicationid);
                            if (affiliations.size() == 0) {
                                continue;
                            }
                            if (std::find(affiliations.begin(), affiliations.end(), NO_AFFILIATION) != affiliations.end())
                            {
                                errorset.insert("get_affiliations() returned error NO_AFFILIATION");
                                continue;
                            }
                            std::vector<Coord> aff_coords = {};
                            aff_coords.reserve(affiliations.size());
                            std::transform(affiliations.begin(),affiliations.end(),std::back_inserter(aff_coords),[this](auto id){return mainprg_.ds_.get_affiliation_coord(id);});
                            if (std::find(aff_coords.begin(),aff_coords.end(),NO_COORD) != aff_coords.end())
                            {
                                errorset.insert("get_affiliations() returned error NO_AFFILIATION");
                                continue;
                            }
                            auto pen = QPen(publicationcolor);
                            pen.setWidth(0); // "Cosmetic" pen
                            Coord centercoord = NO_COORD;
                            long long int x = std::accumulate(aff_coords.begin(), aff_coords.end(), 0, [](auto& coord1, auto& coord2){return coord1 + coord2.x;});
                            x /= aff_coords.size();
                            long long int y = std::accumulate(aff_coords.begin(), aff_coords.end(), 0, [](auto& coord1, auto& coord2){return coord1 + coord2.y;});
                            y /= aff_coords.size();
                            centercoord.x = x;
                            centercoord.y = y;
                            for (auto& coord : aff_coords)
                            {
                                if (centercoord != NO_COORD)
                                {
                                    QLineF line(QPointF(20*centercoord.x, -20*centercoord.y), QPointF(20*coord.x, -20*coord.y));
                                    auto lineitem = gscene_->addLine(line, pen);
                                    lineitem->setFlag(QGraphicsItem::ItemIsSelectable);
                                    lineitem->setData(0, QVariant::fromValue(publicationid));
                                    lineitem->setZValue(publicationzvalue);
                                }
                            }

                        }
                    }
                    catch (NotImplemented const& e)
                    {
                        errorset.insert(std::string("NotImplemented while updating graphics: ") + e.what());
                        std::cerr << std::endl << "NotImplemented while updating graphics: " << e.what() << std::endl;
                    }
                }
            }
            catch (NotImplemented const& e)
            {
                errorset.insert(std::string("NotImplemented while updating graphics: ") + e.what());
                std::cerr << std::endl << "NotImplemented while updating graphics: " << e.what() << std::endl;
            }
        }

        if (ui->connections_checkbox->isChecked()) {
            try {
                auto connections = mainprg_.ds_.get_all_connections();
                QColor connectioncolor = Qt::magenta;
                auto min_weight = ui->min_weight_spinbox->value();
                for (auto& connection : connections){
                    if (connection.weight<min_weight){
                        continue;
                    }
                    auto start_coord = mainprg_.ds_.get_affiliation_coord(connection.aff1);
                    auto end_coord = mainprg_.ds_.get_affiliation_coord(connection.aff2);
                    auto pen = QPen(connectioncolor);
                    pen.setWidth(0); // "Cosmetic" pen
                    QLineF line(QPointF(20*start_coord.x, -20*start_coord.y), QPointF(20*end_coord.x, -20*end_coord.y));
                    auto lineitem = gscene_->addLine(line, pen);
                    lineitem->setZValue(0);
                }
                if(!result_connections.empty()){
                    for (auto& connection: result_connections){
                        auto start_coord = mainprg_.ds_.get_affiliation_coord(connection.aff1);
                        auto end_coord = mainprg_.ds_.get_affiliation_coord(connection.aff2);
                        auto pen = QPen(Qt::green);
                        if (connection.weight<min_weight){
                            pen.setStyle(Qt::DashLine);
                        }
                        pen.setWidth(0); // "Cosmetic" pen
                        QLineF line(QPointF(20*start_coord.x, -20*start_coord.y), QPointF(20*end_coord.x, -20*end_coord.y));
                        auto lineitem = gscene_->addLine(line, pen);
                        lineitem->setZValue(5);
                    }
                }
            } catch (NotImplemented const& e) {
                errorset.insert(std::string("NotImplemented while updating graphics: ") + e.what());
                std::cerr << std::endl << "NotImplemented while updating graphics: " << e.what() << std::endl;
            }
        }

        #ifdef WORLDMAP_HH
        if (ui->world_map_checkbox->isChecked()) {
            auto pen = QPen(Qt::red);
            pen.setWidth(0); // "Cosmetic" pen
            for (auto& polygon : worldmap::POLYGONS) {
                #if QT_VERSION >= QT_VERSION_CHECK(6,0,0)
                QList<QPointF> list_of_points;
                #else
                QVector<QPointF> list_of_points;
                #endif
                std::for_each(polygon.begin(),polygon.end(),[&list_of_points](auto& pair){list_of_points.push_back(QPointF(pair.first,pair.second));});
                QPolygonF qpolygon(list_of_points);
                auto polygonitem = gscene_->addPolygon(qpolygon,pen);
                polygonitem->setZValue(0);
            }
        }
        #endif
    }
    catch (NotImplemented const& e)
    {
        errorset.insert(std::string("NotImplemented while updating graphics: ") + e.what());
        std::cerr << std::endl << "NotImplemented while updating graphics: " << e.what() << std::endl;
    }

    if (!errorset.empty())
    {
        std::ostringstream errorstream;
        for (auto const& errormsg : errorset)
        {
            errorstream << "Error from GUI: " << errormsg << std::endl;
        }
        output_text(errorstream);
        output_text_end();
    }
}

void MainWindow::output_text(ostringstream& output)
{
    string outstr = output.str();
    if (!outstr.empty())
    {
        if (outstr.back() == '\n') { outstr.pop_back(); } // Remove trailing newline
        ui->output->appendPlainText(QString::fromStdString(outstr));
        ui->output->ensureCursorVisible();
        ui->output->repaint();
    }

    output.str(""); // Clear the stream, because it has already been output
}

void MainWindow::output_text_end()
{
    ui->output->moveCursor(QTextCursor::End);
    ui->output->ensureCursorVisible();
    ui->output->repaint();
}

bool MainWindow::check_stop_pressed() const
{
    QCoreApplication::processEvents();
    return stop_pressed_;
}

void MainWindow::execute_line()
{
    auto line = ui->lineEdit->text();
    clear_input_line();
    ui->output->appendPlainText(QString::fromStdString(MainProgram::PROMPT)+line);

    ui->execute_button->setEnabled(false);
    ui->stop_button->setEnabled(true);
    stop_pressed_ = false;

    ostringstream output;
    bool cont = mainprg_.command_parse_line(line.toStdString(), output);
    ui->lineEdit->clear();
    output_text(output);
    output_text_end();

    ui->stop_button->setEnabled(false);
    ui->execute_button->setEnabled(true);
    stop_pressed_ = false;

    ui->lineEdit->setFocus();

    update_view();

    if (!cont)
    {
        close();
    }
}

void MainWindow::cmd_selected(int idx)
{
    ui->lineEdit->insert(QString::fromStdString(mainprg_.cmds_[idx].cmd+" "));
    ui->cmd_info_text->setText(QString::fromStdString(mainprg_.cmds_[idx].cmd+" "+mainprg_.cmds_[idx].info));

    ui->lineEdit->setFocus();
}

void MainWindow::number_selected(QString const& number)
{
    ui->lineEdit->insert(number+" ");

    ui->lineEdit->setFocus();
}

void MainWindow::select_file()
{
    QFileDialog dialog(this, "Select file", QDir::currentPath(), "Command list (*.txt)");
    dialog.setFileMode(QFileDialog::ExistingFiles);
    if (dialog.exec())
    {
        auto filenames = dialog.selectedFiles();
        for (auto& i : filenames)
        {
            auto filename = QDir("").relativeFilePath(i);
            ui->lineEdit->insert("\""+filename+"\" ");
        }
    }

    ui->lineEdit->setFocus();
}

void MainWindow::clear_input_line()
{
    ui->cmd_info_text->clear();
    ui->lineEdit->clear();
    ui->lineEdit->setFocus();
}

void MainWindow::fit_view()
{
    ui->graphics_view->fitInView(gscene_->itemsBoundingRect(), Qt::KeepAspectRatio);
}

void MainWindow::scene_selection_change()
{
    auto items = gscene_->selectedItems();
    if (!items.empty())
    {
        for (auto i : items)
        {
            auto data = i->data(0);
            if (!selection_clear_in_progress)
            {
                ostringstream output;
                output << "*click* ";
                std::string ui_str;
                if (data.canConvert<AffiliationID>())
                {
                    auto affiliationid = data.value<AffiliationID>();
                    ui_str = mainprg_.print_affiliation(affiliationid, output);
                }
                if (data.canConvert<PublicationID>())
                {
                    auto publicationid = data.value<PublicationID>();
                    ui_str = mainprg_.print_publication(publicationid, output);
                }
                else if (data.canConvert<Coord>())
                {
                    auto coord = data.value<Coord>();
                    ui_str = mainprg_.print_coord(coord, output);
                }
                if (!ui_str.empty()) { ui->lineEdit->insert(QString::fromStdString(ui_str+" ")); }
                output_text(output);
                output_text_end();
            }
            i->setSelected(false);
            selection_clear_in_progress = !selection_clear_in_progress;
        }
    }
}

void MainWindow::clear_selection()
{
    gscene_->clearSelection();
}



// Originally in main.cc
#include <QApplication>

int main(int argc, char *argv[])
{
    if (argc == 1)
    {
        QApplication a(argc, argv);
        MainWindow w;
        w.show();

        auto status = a.exec();
        std::cerr << "Program ended normally." << std::endl;
        return status;
    }
    else
    {
        return MainProgram::mainprogram(argc, argv);
    }
}

#else

int main(int argc, char *argv[])
{
    return MainProgram::mainprogram(argc, argv);
}
#endif
