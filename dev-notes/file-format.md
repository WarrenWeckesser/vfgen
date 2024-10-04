VFGEN needs an alternative to the XML file format.

For the current capabilities of VFGEN, the essential information is

vector_field
    name                    (id)
    independent_variable    (id)
    description             (text)

constant
    name                    (id)
    value                   (formula)
    description             (text)
    latex                   (text)

parameter
    name                    (id)
    description             (text)
    default_value           (formula)
    latex                   (text)

expression
    name                    (id)
    formula                 (formula)
    description             (text)
    latex                   (text)

state_variable
    name                    (id)
    derivative              (formula)
    description             (text)
    period_from             (formula)
    period_to               (formula)
    default_ic              (formula)
    default_history         (formula)
    latex                   (text)

function
    name                    (id)
    formula                 (formula)
    description             (text)
