#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""My first try of building a small dialog for editing filenames."""

try:
    from enthought.traits.api import HasTraits, List, Str    
    from enthought.traits.ui.api import View, Item, ListEditor, Group
    from enthought.traits.ui.menu import OKButton, CancelButton, RevertButton, UndoButton
except ImportError:
    from traits.api import HasTraits, List, Str    
    from traitsui.api import View, Item, ListEditor, Group
    from traitsui.menu import OKButton, CancelButton, RevertButton, UndoButton

class ChannelNamesEditor ( HasTraits ):
    """ChannelNamesEditor
    Used to edit the channel-names of EEG files."""
    # THE trait-attribute to change
    channel_names = List( Str )
    view = View(
        Group(
            Item( 'channel_names@',
                  show_label = False,
                  editor = ListEditor(style="text")
            ),
            label="Channel names:"
        ), 
        title     = 'Editing channel-names',
        width     = 0.5,
        height    = 0.6,
        resizable = True,
        buttons=[OKButton, CancelButton, RevertButton, UndoButton]
    )

if __name__ == "__main__":
    example_names = ["Fp1","Fp2","Fpz","F3","Fz","F4"]
    print "Editing names:", example_names
    ce = ChannelNamesEditor(channel_names=example_names)
    ce.configure_traits()
    example_names = ce.channel_names
    print "Now the names are:", example_names
    
    