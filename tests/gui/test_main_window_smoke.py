import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def test_main_window_constructs(qtbot):
    from bvidfe.gui.main_window import BvidMainWindow

    w = BvidMainWindow()
    qtbot.addWidget(w)
    assert w.isVisible() is False  # we didn't show it
    assert w.windowTitle() == "BVID-FE"
