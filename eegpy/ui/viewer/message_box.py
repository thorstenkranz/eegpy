from enthought.traits.api import *
from enthought.traits.ui.api import View, Item
from enthought.traits.ui.menu import NoButtons, OKButton, ToolBar, Action

                
class Message(HasTraits):
    message = Str()
    view = View( Item('message',
                      style='readonly',
                      show_label=False
                 ),
                 buttons=[OKButton],
                 title="Error",
                 width=0.4,
                 kind="modal",
               ) 
    