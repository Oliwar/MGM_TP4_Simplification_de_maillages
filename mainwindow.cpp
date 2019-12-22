#include "mainwindow.h"
#include "ui_mainwindow.h"

/* **** début de la partie à compléter **** */
/**
 * Affiche l'arête actuellement sélectionné
 * @brief MainWindow::showEdgeSelection
 * @param _mesh
 */
void MainWindow::showEdgeSelection(MyMesh* _mesh)
{
    // on réinitialise les couleurs de tout le maillage
    resetAllColorsAndThickness(_mesh);

    /* **** à compléter ! (Partie 1) ****
     * cette fonction utilise la variables de sélection edgeSelection qui est l'ID de
     * l'arête sélectionnée et qui est égale à -1 si la sélection est vide
     */
    if(edgeSelection >= 0 && edgeSelection < _mesh->n_edges()){
        EdgeHandle eh = _mesh->edge_handle(edgeSelection); // EdgeHandle d'indice i
        _mesh->set_color(eh, MyMesh::Color(0, 255, 0));
        _mesh->data(eh).thickness = 4;

        HalfedgeHandle heh = _mesh->halfedge_handle(eh, 0);
        VertexHandle vh1 = _mesh->to_vertex_handle(heh);
        _mesh->set_color(vh1, MyMesh::Color(255, 0, 0));
        _mesh->data(vh1).thickness = 8;

        VertexHandle vh2 = _mesh->from_vertex_handle(heh);
        _mesh->set_color(vh2, MyMesh::Color(255, 0, 0));
        _mesh->data(vh2).thickness = 8;
    }

    // on affiche le nouveau maillage
    displayMesh(_mesh);
}

/**
 * Marque les arêtes à supprimer
 * @brief MainWindow::collapseEdge
 * @param _mesh
 * @param edgeID
 */
void MainWindow::collapseEdge(MyMesh* _mesh, int edgeID)
{
    /* **** à compléter ! (Partie 1) ****
     * cette fonction utilise l'opérateur collapse pour supprimer l'arête d'index edgeID
     * Attention à ne pas oublier garbage_collection() !
     */
    if(edgeID >= 0 && edgeID < _mesh->n_edges()){
        EdgeHandle eh = _mesh->edge_handle(edgeID);

        HalfedgeHandle heh = _mesh->halfedge_handle(eh, 0);
        if(_mesh->is_collapse_ok(heh)){
            VertexHandle vh = _mesh->to_vertex_handle(heh);
            MyMesh::Point middle = _mesh->calc_edge_midpoint(eh);

            _mesh->collapse(heh);
            _mesh->set_point(vh, middle);

            // permet de nettoyer le maillage et de garder la cohérence des indices après un collapse
            //_mesh->garbage_collection();
        }
    }

}

// fonction pratique pour faire des tirages aléatoires
int randInt(int low, int high){return qrand() % ((high + 1) - low) + low;}

/**
 * Suppression aléatoire des arêtes
 * @brief MainWindow::random_delete
 * @param _mesh
 * @param percent
 */
void MainWindow::random_delete(MyMesh* _mesh, int percent){
    unsigned nb_aretes_del = _mesh->n_edges() * percent / 100;

    qDebug() << "nb_totales edges " << _mesh->n_edges();
    qDebug() << "nb aretes a del " << nb_aretes_del;

    unsigned rng;
    EdgeHandle eh;
    HalfedgeHandle heh;
    VertexHandle vh_from;
    VertexHandle vh_to;
    for(unsigned i = 0; i < nb_aretes_del;){
        while(true){
            rng = randInt(0, _mesh->n_edges() - 1);
            eh = _mesh->edge_handle(rng);

            if(_mesh->status(eh).deleted()) continue;

            heh = _mesh->halfedge_handle(eh, 0);
            if(!_mesh->is_collapse_ok(heh)) continue;

            vh_from = _mesh->from_vertex_handle(heh);
            vh_to = _mesh->to_vertex_handle(heh);
            if(!_mesh->is_boundary(eh))
                if(_mesh->is_boundary(vh_from) || _mesh->is_boundary(vh_to))
                    continue;

            collapseEdge(_mesh, eh.idx());

            if(!_mesh->is_boundary(eh)) i += 3;
            else i += 2;

            break;
        }
    }
    _mesh->garbage_collection();
    qDebug() << "nb aretes finale " << _mesh->n_edges();
}

/**
 * Suppression des arêtes de plus petite taille en premier
 * @brief MainWindow::size_delete
 * @param _mesh
 * @param percent
 */
void MainWindow::size_delete(MyMesh* _mesh, int percent){
    unsigned nb_aretes_del = _mesh->n_edges() * percent / 100;

    qDebug() << "nb_totales edges " << _mesh->n_edges();
    qDebug() << "nb aretes a del " << nb_aretes_del;

    std::vector<float> size = std::vector<float>(_mesh->n_edges(), 0.0f);

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
        size.at(curEdge->idx()) = _mesh->calc_edge_length(*curEdge);

    unsigned minElementIndex;
    VertexHandle vh_to;
    VertexHandle vh_from;
    HalfedgeHandle heh;
    EdgeHandle eh;
    for(size_t i = 0; i < nb_aretes_del;){
        minElementIndex = std::min_element(size.begin(), size.end()) - size.begin();

        eh = _mesh->edge_handle(minElementIndex);
        size.at(eh.idx()) = 1000.0f;

        heh = _mesh->halfedge_handle(eh, 0);
        if(!_mesh->is_collapse_ok(heh)) continue;

        vh_from = _mesh->from_vertex_handle(heh);
        vh_to = _mesh->to_vertex_handle(heh);
        if(!_mesh->is_boundary(eh))
            if(_mesh->is_boundary(vh_from) || _mesh->is_boundary(vh_to)) continue;

        collapseEdge(_mesh, eh.idx());

        for (MyMesh::VertexEdgeIter curEdge = _mesh->ve_iter(vh_to); curEdge.is_valid(); curEdge ++)        
            size.at(curEdge->idx()) = _mesh->calc_edge_length(*curEdge);

        if(!_mesh->is_boundary(eh)) i += 3;
        else i += 2;
    }
    _mesh->garbage_collection();
    qDebug() << "nb aretes finale " << _mesh->n_edges();
}

/**
 * Suppression des arêtes de plus grand angle en premier
 * @brief MainWindow::angle_delete
 * @param _mesh
 * @param percent
 */
void MainWindow::angle_delete(MyMesh* _mesh, int percent){
    unsigned nb_aretes_del = static_cast<unsigned>(_mesh->n_edges() * static_cast<unsigned>(percent) / 100);

    qDebug() << "Nombre totale d'aretes de depart : " << _mesh->n_edges();
    qDebug() << "Nombre d'aretes à supprimer : " << nb_aretes_del;

    std::vector<float> angle = std::vector<float>(_mesh->n_edges(), 0.0f);

    float angle_dihedral;
    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++){
        if(!_mesh->is_boundary(*curEdge)){
            angle_dihedral = _mesh->calc_dihedral_angle(*curEdge);
            if(angle_dihedral == 0.0f) angle_dihedral = static_cast<float>(M_PI);
            angle.at(static_cast<unsigned long>(curEdge->idx())) = std::abs(angle_dihedral);
        }
    }

    unsigned maxElementIndex;
    VertexHandle vh_to;
    VertexHandle vh_from;
    HalfedgeHandle heh;
    EdgeHandle eh;
    for(unsigned i = 0; i < nb_aretes_del;){
        maxElementIndex = static_cast<unsigned>(std::max_element(angle.begin(), angle.end()) - angle.begin());

        eh = _mesh->edge_handle(maxElementIndex);
        angle.at(static_cast<unsigned long>(eh.idx())) = -1.0f;

        heh = _mesh->halfedge_handle(eh, 0);
        if(!_mesh->is_collapse_ok(heh)) continue;

        vh_from = _mesh->from_vertex_handle(heh);
        vh_to = _mesh->to_vertex_handle(heh);
        if(!_mesh->is_boundary(eh))
            if(_mesh->is_boundary(vh_from) || _mesh->is_boundary(vh_to)) continue;

        collapseEdge(_mesh, eh.idx());

        if(!_mesh->is_boundary(eh)) i += 3;
        else i += 2;

        /*for (MyMesh::VertexOHalfedgeIter curHalfEdge = _mesh->voh_iter(vh_to); curHalfEdge.is_valid(); curHalfEdge ++){
            eh = _mesh->edge_handle(*curHalfEdge);
            if(!_mesh->is_boundary(eh)){
                angle_dihedral = _mesh->calc_dihedral_angle(eh);
                if(angle_dihedral == 0.0f) angle_dihedral = static_cast<float>(M_PI);
                angle.at(static_cast<unsigned long>(eh.idx())) = std::abs(angle_dihedral);
            } else angle.at(static_cast<unsigned long>(eh.idx())) = 0.0f;

            eh = _mesh->edge_handle(_mesh->next_halfedge_handle(*curHalfEdge));
            if(!_mesh->is_boundary(eh)){
                angle_dihedral = _mesh->calc_dihedral_angle(eh);
                if(angle_dihedral == 0.0f) angle_dihedral = static_cast<float>(M_PI);
                angle.at(static_cast<unsigned long>(eh.idx())) = std::abs(angle_dihedral);
            } else angle.at(static_cast<unsigned long>(eh.idx())) = 0.0f;
        }*/
    }
    _mesh->garbage_collection();
    qDebug() << "Nombre d'aretes finales : " << _mesh->n_edges();
}

/**
 * Suppression des arêtes avec une planeite la plus grande en premier
 * @brief MainWindow::flatness_delete
 * @param _mesh
 * @param percent
 */
void MainWindow::flatness_delete(MyMesh* _mesh, int percent){
    unsigned nb_aretes_del = static_cast<unsigned>(_mesh->n_edges() * static_cast<unsigned>(percent) / 100);

    qDebug() << "Nombre totale d'aretes de depart : " << _mesh->n_edges();
    qDebug() << "Nombre d'aretes à supprimer : " << nb_aretes_del;

    std::vector<float> planeite = std::vector<float>(_mesh->n_edges(), 0.0f);

    float angle_dihedral;
    VertexHandle vh_to;
    VertexHandle vh_from;
    HalfedgeHandle heh;
    float cpt;
    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++){
        heh = _mesh->halfedge_handle(*curEdge, 0);
        vh_to = _mesh->to_vertex_handle(heh);
        vh_from = _mesh->from_vertex_handle(heh);

        cpt = 0;
        angle_dihedral = 0;
        for (MyMesh::VertexEdgeIter curVEdge = _mesh->ve_iter(vh_to); curVEdge.is_valid(); curVEdge ++){
            if(!_mesh->is_boundary(*curVEdge)){
                if(_mesh->calc_dihedral_angle(*curVEdge) == 0.0f) angle_dihedral += static_cast<float>(M_PI);
                else angle_dihedral += _mesh->calc_dihedral_angle(*curVEdge);
                cpt++;
            }
        }
        for (MyMesh::VertexEdgeIter curVEdge = _mesh->ve_iter(vh_from); curVEdge.is_valid(); curVEdge ++){
            if(!_mesh->is_boundary(*curVEdge)){
                if(_mesh->calc_dihedral_angle(*curVEdge) == 0.0f) angle_dihedral += static_cast<float>(M_PI);
                else angle_dihedral += _mesh->calc_dihedral_angle(*curVEdge);
                cpt++;
            }
        }
        planeite.at(static_cast<unsigned long>(curEdge->idx())) = std::abs(angle_dihedral)/cpt;
    }

    unsigned maxElementIndex;
    EdgeHandle eh;
    for(unsigned i = 0; i < nb_aretes_del;){
        maxElementIndex = static_cast<unsigned>(std::max_element(planeite.begin(), planeite.end()) - planeite.begin());

        eh = _mesh->edge_handle(maxElementIndex);
        planeite.at(static_cast<unsigned long>(eh.idx())) = -1.0f;

        heh = _mesh->halfedge_handle(eh, 0);
        if(!_mesh->is_collapse_ok(heh)) continue;

        //Garde l'integrité des bordures
        vh_to = _mesh->to_vertex_handle(heh);
        vh_from = _mesh->from_vertex_handle(heh);
        if(!_mesh->is_boundary(eh))
            if(_mesh->is_boundary(vh_from) || _mesh->is_boundary(vh_to)) continue;

        collapseEdge(_mesh, eh.idx());

        /*VertexHandle vh_to2;
        for (MyMesh::VertexVertexIter curVertex = _mesh->vv_iter(vh_to); curVertex.is_valid(); curVertex ++)
        {
            for (MyMesh::VertexEdgeIter curEdge = _mesh->ve_iter(*curVertex); curEdge.is_valid(); curEdge ++){
                heh = _mesh->halfedge_handle(*curEdge, 0);
                vh_to2 = _mesh->to_vertex_handle(heh);
                vh_from = _mesh->from_vertex_handle(heh);

                cpt = 0;
                angle_dihedral = 0;
                for (MyMesh::VertexEdgeIter curVEdge = _mesh->ve_iter(vh_to2); curVEdge.is_valid(); curVEdge ++){
                    if(!_mesh->is_boundary(*curVEdge)){
                        if(_mesh->calc_dihedral_angle(*curVEdge) == 0.0f) angle_dihedral += static_cast<float>(M_PI);
                        else angle_dihedral += _mesh->calc_dihedral_angle(*curVEdge);
                        cpt++;
                    }
                }
                for (MyMesh::VertexEdgeIter curVEdge = _mesh->ve_iter(vh_from); curVEdge.is_valid(); curVEdge ++){
                    if(!_mesh->is_boundary(*curVEdge)){
                        if(_mesh->calc_dihedral_angle(*curVEdge) == 0.0f) angle_dihedral += static_cast<float>(M_PI);
                        else angle_dihedral += _mesh->calc_dihedral_angle(*curVEdge);
                        cpt++;
                    }
                }
                planeite.at(static_cast<unsigned long>(curEdge->idx())) = std::abs(angle_dihedral)/cpt;
            }
        }*/

        if(!_mesh->is_boundary(eh)) i += 3;
        else i += 2;
    }
    _mesh->garbage_collection();
    qDebug() << "nb aretes finale " << _mesh->n_edges();
}

void MainWindow::decimation(MyMesh* _mesh, int percent, QString method)
{
    /* **** à compléter ! (Partie 2 et 3) ****
     * Cette fonction supprime des arêtes jusqu'à atteindre un pourcentage d'arêtes restantes, selon un critère donné
     * percent : pourcentage de l'objet à garder
     * method  : la méthode à utiliser parmis : "Aléatoire", "Par taille", "Par angle", "Par planéité"
     */

    if(method == "Aléatoire")
    {
        random_delete(_mesh, percent);
    }
    else if(method == "Par taille")
    {
        /*unsigned nb_aretes_del = _mesh->n_edges() * percent / 100;

        qDebug() << "nb_totales edges " << _mesh->n_edges();
        qDebug() << "nb aretes a del " << nb_aretes_del;

        std::vector<Edges> *size = new std::vector<Edges>(_mesh->n_edges());
        for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++){
            Edges edges;
            edges.id = curEdge->idx();
            edges.size = _mesh->calc_edge_length(*curEdge);
            size->push_back(edges);
        }

        std::sort(size->begin(), size->end());

        VertexHandle vh;
        HalfedgeHandle heh;
        EdgeHandle eh;
        qDebug() << "init ok";
        for(size_t i = 0; i < nb_aretes_del;){
            eh = _mesh->edge_handle(size->at(0).id);
            heh = _mesh->halfedge_handle(eh, 0);
            if(!_mesh->is_collapse_ok(heh)){
                size->erase(size->begin());
                continue;
            }

            vh = _mesh->to_vertex_handle(heh);

            collapseEdge(_mesh, eh.idx());

            qDebug() << "avant";
            for (MyMesh::VertexEdgeIter curEdge = _mesh->ve_iter(vh); curEdge.is_valid(); curEdge ++)
            {
                if(_mesh->status(*curEdge).deleted()){
                    qDebug() << "else";
                    for(size_t i = 0; i < size->size(); ++i){
                        qDebug() << "if";
                        if(size->at(i).id == curEdge->idx()){
                            size->erase(size->begin() + i);
                            break;
                        }
                    }
                    qDebug() << "fin else";
                }
                else{
                    for(size_t i = 0; i < size->size(); i++){
                        qDebug() << "if";
                        if(size->at(i).id == curEdge->idx()){
                            size->at(curEdge->idx()).size = _mesh->calc_edge_length(*curEdge);
                            break;
                        }
                    }
                }
            }
            qDebug() << "apres";

            if(!_mesh->is_boundary(eh)) i += 3;
            else i += 2;
        }
        _mesh->garbage_collection();
        qDebug() << "nb aretes finale " << _mesh->n_edges();
*/
        size_delete(_mesh, percent);
    }
    else if(method == "Par angle")
    {
        angle_delete(_mesh, percent);
    }
    else if(method == "Par planéité")
    {
        flatness_delete(_mesh, percent);
    }
    else
    {
        qDebug() << "Méthode inconnue !!!";
    }

}

/* **** début de la partie boutons et IHM **** */
void MainWindow::updateEdgeSelectionIHM()
{
    /* **** à compléter ! (Partie 3) ****
     * Cette fonction met à jour l'interface, les critères pourrons être affichés dans la zone de texte pour les vérifier
     */

    QString infos = "";
    infos = infos + "Surface : " + QString::number(0) + "\n";
    infos = infos + "C1 : " + QString::number(0) + "\n";
    infos = infos + "C2 : " + QString::number(0) + "\n";
    infos = infos + "C3 : " + QString::number(0) + "\n";
    ui->infoEdgeSelection->setPlainText(infos);

    ui->labelEdge->setText(QString::number(edgeSelection));

    // on montre la nouvelle sélection
    showEdgeSelection(&mesh);
}
/* **** fin de la partie à compléter **** */

void MainWindow::on_pushButton_edgeMoins_clicked()
{
    // mise à jour de l'interface
    edgeSelection = edgeSelection - 1;
    updateEdgeSelectionIHM();
}

void MainWindow::on_pushButton_edgePlus_clicked()
{
    // mise à jour de l'interface
    edgeSelection = edgeSelection + 1;
    updateEdgeSelectionIHM();
}

void MainWindow::on_pushButton_delSelEdge_clicked()
{
    // on supprime l'arête d'indice edgeSelection
    collapseEdge(&mesh, edgeSelection);

    // permet de nettoyer le maillage et de garder la cohérence des indices après un collapse
    mesh.garbage_collection();

    // on actualise la sélection
    showEdgeSelection(&mesh);
}

void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);
}

void MainWindow::on_pushButton_decimate_clicked()
{
    decimation(&mesh, ui->horizontalSlider->value(), ui->comboBox->currentText());
    displayMesh(&mesh);
}
/* **** fin de la partie boutons et IHM **** */



/* **** fonctions supplémentaires **** */
// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }
}

// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh* _mesh, bool isTemperatureMap, float mapRange)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    int i = 0;

    if(isTemperatureMap)
    {
        QVector<float> values;

        if(mapRange == -1)
        {
            for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
                values.append(fabs(_mesh->data(*curVert).value));
            qSort(values);
            mapRange = values.at(values.size()*0.8);
            qDebug() << "mapRange" << mapRange;
        }

        float range = mapRange;
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;

        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }
    else
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }


    ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;

    GLuint* linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat* linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat* linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    i = 0;
    QHash<float, QList<int> > edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if(t > 0)
        {
            if(!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int> > it(edgesIDbyThickness);
    QList<QPair<float, int> > edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for(int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh1)[0];
            linesVerts[3*i+1] = _mesh->point(vh1)[1];
            linesVerts[3*i+2] = _mesh->point(vh1)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh2)[0];
            linesVerts[3*i+1] = _mesh->point(vh2)[1];
            linesVerts[3*i+2] = _mesh->point(vh2)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint* pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat* pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat* pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int> > vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if(t > 0)
        {
            if(!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }
    QHashIterator<float, QList<int> > vitt(vertsIDbyThickness);
    QList<QPair<float, int> > vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for(int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3*i+0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3*i+1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3*i+2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3*i+0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3*i+1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3*i+2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    vertexSelection = -1;
    edgeSelection = -1;
    faceSelection = -1;

    modevoisinage = false;

    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

