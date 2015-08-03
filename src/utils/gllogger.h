#ifndef GLLOGGER_H
#define GLLOGGER_H

#include <QOpenGLDebugLogger>


namespace utils {

namespace ogl {

class Logger : public QOpenGLDebugLogger
{
    Q_OBJECT

  public:
    explicit Logger(QObject *parent = 0)
      : QOpenGLDebugLogger(parent)
    { }

    bool init(void);

  private slots:
    void handleLoggedMessage(const QOpenGLDebugMessage & msg);
};

} // End of namespace ogl

} // End of namespace utils

#endif // GLLOGGER_H
