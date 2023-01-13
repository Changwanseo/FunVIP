# ete3 patch
from __future__ import absolute_import
from __future__ import print_function
import ete3

# Patch for ete3 3.1.2 colliding with higher version PyQt5
def save(scene, imgName, w=None, h=None, dpi=90, take_region=False, units="px"):

    import colorsys
    import random
    import re
    import types
    from sys import stderr

    from PyQt5 import QtGui, QtCore
    from PyQt5.QtCore import (
        Qt,
        QPointF,
        QRect,
        QRectF,
        QBuffer,
        QByteArray,
        QThread,
        QIODevice,
        QMetaObject,
        QModelIndex,
        QObject,
        QRegExp,
        QSize,
        QSizeF,
        QVariant,
    )  # QString
    from PyQt5.QtSvg import QGraphicsSvgItem, QSvgGenerator
    from PyQt5.QtOpenGL import QGLFormat, QGLWidget
    from PyQt5.QtPrintSupport import QPrinter
    from PyQt5.QtWidgets import (
        QAction,
        QApplication,
        QCheckBox,
        QWidget,
        QColorDialog,
        QComboBox,
        QDialog,
        QDialogButtonBox,
        QFileDialog,
        QGraphicsEllipseItem,
        QGraphicsItem,
        QGraphicsItemGroup,
        QGraphicsLineItem,
        QGraphicsPathItem,
        QGraphicsPixmapItem,
        QGraphicsPolygonItem,
        QGraphicsRectItem,
        QGraphicsScene,
        QGraphicsSimpleTextItem,
        QGraphicsTextItem,
        QGraphicsView,
        QInputDialog,
        QItemDelegate,
        QLabel,
        QLineEdit,
        QListWidget,
        QMainWindow,
        QMenu,
        QMenuBar,
        QMessageBox,
        QPushButton,
        QSplitter,
        QStatusBar,
        QTableView,
        QTextEdit,
        QToolBar,
        QVBoxLayout,
        QWidget,
    )

    from PyQt5.QtGui import (
        QBrush,
        QColor,
        QCursor,
        QFont,
        QFontMetrics,
        QIcon,
        QImage,
        QPainter,
        QPainterPath,
        QPen,
        QPixmap,
        QPolygonF,
        QRadialGradient,
        QRegExpValidator,
        QStandardItemModel,
        QTransform,
    )

    QtCore.pyqtSignature = QtCore.pyqtSlot

    from ete3.treeview.svg_colors import SVG_COLORS, COLOR_SCHEMES

    import time

    ipython_inline = False
    if imgName == "%%inline":
        ipython_inline = True
        ext = "PNG"
    elif imgName == "%%inlineSVG":
        ipython_inline = True
        ext = "SVG"
    elif imgName.startswith("%%return"):
        try:
            ext = imgName.split(".")[1].upper()
        except IndexError:
            ext = "SVG"
        imgName = "%%return"
    else:
        ext = imgName.split(".")[-1].upper()

    main_rect = scene.sceneRect()
    aspect_ratio = main_rect.height() / main_rect.width()

    # auto adjust size
    if not w and not h:
        units = "px"
        w = main_rect.width()
        h = main_rect.height()
        ratio_mode = Qt.KeepAspectRatio
    elif w and h:
        ratio_mode = Qt.IgnoreAspectRatio
    elif h is None:
        h = w * aspect_ratio
        ratio_mode = Qt.KeepAspectRatio
    elif w is None:
        w = h / aspect_ratio
        ratio_mode = Qt.KeepAspectRatio

    # Adjust to resolution
    if units == "mm":
        if w:
            w = w * 0.0393700787 * dpi
        if h:
            h = h * 0.0393700787 * dpi
    elif units == "in":
        if w:
            w = w * dpi
        if h:
            h = h * dpi
    elif units == "px":
        pass
    else:
        raise Exception("wrong unit format")

    x_scale, y_scale = w / main_rect.width(), h / main_rect.height()

    if ext == "SVG":
        svg = QSvgGenerator()
        targetRect = QRectF(0, 0, w, h)
        svg.setSize(QSize(int(w), int(h)))
        svg.setViewBox(targetRect)
        svg.setTitle("Generated with ETE http://etetoolkit.org")
        svg.setDescription("Generated with ETE http://etetoolkit.org")

        if imgName == "%%return":
            ba = QByteArray()
            buf = QBuffer(ba)
            buf.open(QIODevice.WriteOnly)
            svg.setOutputDevice(buf)
        else:
            svg.setFileName(imgName)

        pp = QPainter()
        pp.begin(svg)
        scene.render(pp, targetRect, scene.sceneRect(), ratio_mode)
        pp.end()
        if imgName == "%%return":
            compatible_code = str(ba)
            print("from memory")
        else:
            compatible_code = open(imgName).read()
        # Fix a very annoying problem with Radial gradients in
        # inkscape and browsers...
        compatible_code = compatible_code.replace("xml:id=", "id=")
        compatible_code = re.sub(
            'font-size="(\d+)"', 'font-size="\\1pt"', compatible_code
        )
        compatible_code = compatible_code.replace("\n", " ")
        compatible_code = re.sub("<g [^>]+>\s*</g>", "", compatible_code)
        # End of fix
        if ipython_inline:
            from IPython.core.display import SVG

            return SVG(compatible_code)

        elif imgName == "%%return":
            return x_scale, y_scale, compatible_code
        else:
            open(imgName, "w").write(compatible_code)

    elif ext == "PDF" or ext == "PS":
        if ext == "PS":
            format = QPrinter.PostScriptFormat
        else:
            format = QPrinter.PdfFormat

        printer = QPrinter(QPrinter.HighResolution)
        printer.setResolution(dpi)
        printer.setOutputFormat(format)
        printer.setPageSize(QPrinter.A4)
        printer.setPaperSize(QSizeF(w, h), QPrinter.DevicePixel)
        printer.setPageMargins(0, 0, 0, 0, QPrinter.DevicePixel)

        # pageTopLeft = printer.pageRect().topLeft()
        # paperTopLeft = printer.paperRect().topLeft()
        # For PS -> problems with margins
        # print paperTopLeft.x(), paperTopLeft.y()
        # print pageTopLeft.x(), pageTopLeft.y()
        # print  printer.paperRect().height(),  printer.pageRect().height()
        # topleft =  pageTopLeft - paperTopLeft

        printer.setFullPage(True)
        printer.setOutputFileName(imgName)
        pp = QPainter(printer)
        targetRect = QRectF(0, 0, w, h)
        scene.render(pp, targetRect, scene.sceneRect(), ratio_mode)
    else:
        targetRect = QRectF(0, 0, w, h)
        ii = QImage(w, h, QImage.Format_ARGB32)
        ii.fill(QColor(Qt.white).rgb())
        ii.setDotsPerMeterX(dpi / 0.0254)  # Convert inches to meters
        ii.setDotsPerMeterY(dpi / 0.0254)
        pp = QPainter(ii)
        pp.setRenderHint(QPainter.Antialiasing)
        pp.setRenderHint(QPainter.TextAntialiasing)
        pp.setRenderHint(QPainter.SmoothPixmapTransform)

        scene.render(pp, targetRect, scene.sceneRect(), ratio_mode)
        pp.end()
        if ipython_inline:
            ba = QByteArray()
            buf = QBuffer(ba)
            buf.open(QIODevice.WriteOnly)
            ii.save(buf, "PNG")
            from IPython.core.display import Image

            return Image(ba.data())
        elif imgName == "%%return":
            ba = QByteArray()
            buf = QBuffer(ba)
            buf.open(QIODevice.WriteOnly)
            ii.save(buf, "PNG")
            return x_scale, y_scale, ba.toBase64()
        else:
            ii.save(imgName)

    return w / main_rect.width(), h / main_rect.height()


def patch():
    ete3.treeview.main.save = save
